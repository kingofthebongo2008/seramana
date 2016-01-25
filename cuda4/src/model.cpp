#include <fem.hpp> // Fortran EMulation library of fable module

#include "constants.h"
#include "arrays.h"


namespace hw {


typedef array_1d_mn_cpu<float, constants::n > data_arr_float;

data_arr_float costhe;
data_arr_float sinthe;


using namespace fem::major_types;

struct common_par
{
  int naca;
  float tau;
  float epsmax;
  float ptmax;

  common_par() :
    naca(fem::int0),
    tau(fem::float0),
    epsmax(fem::float0),
    ptmax(fem::float0)
  {}
};

struct common_bod
{
  int nlower;
  int nupper;
  int nodtot;
  arr<float> x;
  arr<float> y;
  common_bod() :
    nlower(fem::int0),
    nupper(fem::int0),
    nodtot(fem::int0),
    x(dimension(100), fem::fill0),
    y(dimension(100), fem::fill0)
  {}
};


struct common_cpd
{
  arr<float> cp;

  common_cpd() :
    cp(dimension(100), fem::fill0)
  {}
};

struct common_commonymous
{
  arr<float> t_lift;

  common_commonymous() :
    t_lift(dimension(10000), fem::fill0)
  {}
};

struct common :
  fem::common,
  common_par,
  common_bod,
  common_cpd,
  common_commonymous
{
  fem::variant_core common_cof;
  fem::cmn_sve cofish_sve;
  fem::cmn_sve veldis_sve;
  fem::cmn_sve gauss_sve;

  common(
    int argc,
    char const* argv[])
  :
    fem::common(argc, argv)
  {}
};

//C
//C**********************************************************
//C
void
naca45(
  common& cmn,
  float const& z,
  float& thick,
  float& camber,
  float& beta)
{
  float& epsmax = cmn.epsmax;
  float& ptmax = cmn.ptmax;
  //
  float dcamdx = fem::float0;
  float w = fem::float0;
  //C
  //C  evaluate thickness and camber
  //C  for naca 4- or 5-digit airfoil
  //C
  thick = 0.0f;
  if (z < 1.e-10f) {
    goto statement_100;
  }
  thick = 5.f * cmn.tau * (.2969f * fem::sqrt(z) - z * (.126f + z * (
    .3537f - z * (.2843f - z * .1015f))));
  statement_100:
  if (epsmax == 0.f) {
    goto statement_130;
  }
  if (cmn.naca > 9999) {
    goto statement_140;
  }
  if (z > ptmax) {
    goto statement_110;
  }
  camber = epsmax / ptmax / ptmax * (2.f * ptmax - z) * z;
  dcamdx = 2.f * epsmax / ptmax / ptmax * (ptmax - z);
  goto statement_120;
  statement_110:
  camber = epsmax / fem::pow2((1.f - ptmax)) * (1.f + z - 2.f *
    ptmax) * (1.f - z);
  dcamdx = 2.f * epsmax / fem::pow2((1.f - ptmax)) * (ptmax - z);
  statement_120:
  beta = fem::atan(dcamdx);
  //C
  return;
  //C
  statement_130:
  camber = 0.0f;
  beta = 0.0f;
  //C
  return;
  //C
  statement_140:
  if (z > ptmax) {
    goto statement_150;
  }
  w = z / ptmax;
  camber = epsmax * w * ((w - 3.f) * w + 3.f - ptmax);
  dcamdx = epsmax * 3.f * w * (1.f - w) / ptmax;
  goto statement_120;
  statement_150:
  camber = epsmax * (1.f - z);
  dcamdx = -epsmax;
  goto statement_120;
  //C
}

//C
//C**********************************************************
//C
void
body(
  common& cmn,
  float& z,
  float const& sign,
  float& x,
  float& y)
{
  //C
  //C  return coordinates of point on body surface
  //C
  //C     z = node spacing parameter
  //C     x,y = cartesian coordinates
  //C     sign = +1. for upper surface, -1. for lower surface
  //C
  if (sign < 0.f) {
    z = 1.f - z;
  }
  float thick = fem::float0;
  float camber = fem::float0;
  float beta = fem::float0;
  naca45(cmn, z, thick, camber, beta);
  x = z - sign * thick * fem::sin(beta);
  y = camber + sign * thick * fem::cos(beta);
  //C
}

//C
//C**********************************************************
//C
void
setup(
  common& cmn)
{
  common_write write(cmn);
  // COMMON bod
  int& nlower = cmn.nlower;
  int& nupper = cmn.nupper;
  int& nodtot = cmn.nodtot;
  arr_ref<float> x(cmn.x, dimension(100));
  arr_ref<float> y(cmn.y, dimension(100));
//  arr_ref<float> costhe(cmn.costhe, dimension(100));
//  arr_ref<float> sinthe(cmn.sinthe, dimension(100));

  const float pi = constant_functions::pi();
  
  //C
  //C  set coordinates of nodes on body surface
  //C
  write(6, "(/,/,/,' body shape',/,/,4x,'x',9x,'y',/)");
  int npoints = nlower;
  float sign = -1.0f;
  int nstart = 0;
  int nsurf = fem::int0;
  int n = fem::int0;
  float fract = fem::float0;
  float z = fem::float0;
  int i = fem::int0;
  FEM_DO_SAFE(nsurf, 1, 2) {
    FEM_DO_SAFE(n, 1, npoints) {
      fract = fem::ffloat(n - 1) / fem::ffloat(npoints);
      z = .5f * (1.f - fem::cos(pi * fract));
      i = nstart + n;
      body(cmn, z, sign, x(i), y(i));
      write(6, "(f8.4,f10.4)"), x(i), y(i);
    }
    npoints = nupper;
    sign = 1.0f;
    nstart = nlower;
  }
  nodtot = nupper + nlower;
  x(nodtot + 1) = x(1);
  y(nodtot + 1) = y(1);
  //C
  //C  set slopes of panels
  //C
  float dx = fem::float0;
  float dy = fem::float0;
  float dist = fem::float0;
  FEM_DO_SAFE(i, 1, nodtot) {
    dx = x(i + 1) - x(i);
    dy = y(i + 1) - y(i);
    dist = fem::sqrt(dx * dx + dy * dy);
    sinthe(i) = dy / dist;
    costhe(i) = dx / dist;
  }
  //C
}

struct cofish_save
{
  fem::variant_bindings cof_bindings;
};

//C
//C**********************************************************
//C
void
cofish(
  common& cmn,
  float const& sinalf,
  float const& cosalf)
{
  FEM_CMN_SVE(cofish);
  int nodtot = cmn.nodtot;
  arr_cref<float> x(cmn.x, dimension(100));
  arr_cref<float> y(cmn.y, dimension(100));

  const float  pi2inv = constant_functions::pi2inv();
  //
  common_variant cof(cmn.common_cof, sve.cof_bindings);
  if (is_called_first_time) {
    using fem::mbr; // member of variant common or equivalence
    {
      mbr<float> a(dimension(101, 111));
      mbr<int> kutta;
      cof.allocate(), a, kutta;
    }
  }
  arr_ref<float, 2> a(cof.bind<float>(), dimension(101, 111));
  int& kutta = cof.bind<int>();
  int j = fem::int0;
  int i = fem::int0;
  float xmid = fem::float0;
  float ymid = fem::float0;
  float flog = fem::float0;
  float ftan = fem::float0;
  float dxj = fem::float0;
  float dxjp = fem::float0;
  float dyj = fem::float0;
  float dyjp = fem::float0;
  float ctimtj = fem::float0;
  float stimtj = fem::float0;
  float b = fem::float0;
  //C
  //C  set coefficients of linear system
  //C
  kutta = nodtot + 1;
  //C
  //C  initialize coefficients
  //C
  FEM_DO_SAFE(j, 1, kutta) {
    a(kutta, j) = 0.0f;
  }
  //C
  //C  set vn=0. at midpoint of ith panel
  //C
  FEM_DO_SAFE(i, 1, nodtot) {
    xmid = .5f * (x(i) + x(i + 1));
    ymid = .5f * (y(i) + y(i + 1));
    a(i, kutta) = 0.0f;
    //C
    //C  find contribution of jth panel
    //C
    FEM_DO_SAFE(j, 1, nodtot) {
      flog = 0.0f;
      ftan = constant_functions::pi();
      if (j == i) {
        goto statement_100;
      }
      dxj = xmid - x(j);
      dxjp = xmid - x(j + 1);
      dyj = ymid - y(j);
      dyjp = ymid - y(j + 1);
      flog = .5f * fem::alog((dxjp * dxjp + dyjp * dyjp) / (dxj *
        dxj + dyj * dyj));
      ftan = fem::atan2(dyjp * dxj - dxjp * dyj, dxjp * dxj + dyjp * dyj);
      statement_100:
      ctimtj = costhe(i) * costhe(j) + sinthe(i) * sinthe(j);
      stimtj = sinthe(i) * costhe(j) - sinthe(j) * costhe(i);
      a(i, j) = pi2inv * (ftan * ctimtj + flog * stimtj);
      b = pi2inv * (flog * ctimtj - ftan * stimtj);
      a(i, kutta) += b;
      if ((i > 1) && (i < nodtot)) {
        goto statement_110;
      }
      //C
      //C  if ith panel touches trailing edge, add contribution
      //C    to kutta condition
      //C
      a(kutta, j) = a(kutta, j) - b;
      a(kutta, kutta) += a(i, j);
      statement_110:;
    }
    //C
    //C  fill in known sides
    //C
    a(i, kutta + 1) = sinthe(i) * cosalf - costhe(i) * sinalf;
  }
  a(kutta, kutta + 1) = -(costhe(1) + costhe(nodtot)) * cosalf - (
    sinthe(1) + sinthe(nodtot)) * sinalf;
  //C
}

struct veldis_save
{
  fem::variant_bindings cof_bindings;
};

//C
//C**********************************************************
//C
void
veldis(
  common& cmn,
  float const& sinalf,
  float const& cosalf)
{
  FEM_CMN_SVE(veldis);
  int& nodtot = cmn.nodtot;
  arr_cref<float> x(cmn.x, dimension(100));
  arr_cref<float> y(cmn.y, dimension(100));

  arr_ref<float> cp(cmn.cp, dimension(100));
  const float pi2inv = constant_functions::pi2inv();
  //
  common_variant cof(cmn.common_cof, sve.cof_bindings);
  if (is_called_first_time) {
    using fem::mbr; // member of variant common or equivalence
    {
      mbr<float> a(dimension(101, 111));
      mbr<int> kutta;
      cof.allocate(), a, kutta;
    }
  }
  arr_cref<float, 2> a(cof.bind<float>(), dimension(101, 111));
  int const& kutta = cof.bind<int>();
  int i = fem::int0;
  arr_1d<150, float> q(fem::fill0);
  float gamma = fem::float0;
  float xmid = fem::float0;
  float ymid = fem::float0;
  float vtang = fem::float0;
  int j = fem::int0;
  float flog = fem::float0;
  float ftan = fem::float0;
  float dxj = fem::float0;
  float dxjp = fem::float0;
  float dyj = fem::float0;
  float dyjp = fem::float0;
  float ctimtj = fem::float0;
  float stimtj = fem::float0;
  float aa = fem::float0;
  float b = fem::float0;
  //C
  //C  compute and print out pressure distribution
  //C
  //C      write(6,1000)
  //C
  //C  retrieve solution from a-matrix
  //C
  FEM_DO_SAFE(i, 1, nodtot) {
    q(i) = a(i, kutta + 1);
  }
  gamma = a(kutta, kutta + 1);
  //C
  //C  find vt and cp at midpoint of ith panel
  //C
  FEM_DO_SAFE(i, 1, nodtot) {
    xmid = .5f * (x(i) + x(i + 1));
    ymid = .5f * (y(i) + y(i + 1));
    vtang = cosalf * costhe(i) + sinalf * sinthe(i);
    //C
    //C  add contributions of jth panel
    //C
    FEM_DO_SAFE(j, 1, nodtot) {
      flog = 0.0f;
      ftan = constant_functions::pi();
      if (j == i) {
        goto statement_100;
      }
      dxj = xmid - x(j);
      dxjp = xmid - x(j + 1);
      dyj = ymid - y(j);
      dyjp = ymid - y(j + 1);
      flog = .5f * fem::alog((dxjp * dxjp + dyjp * dyjp) / (dxj *
        dxj + dyj * dyj));
      ftan = fem::atan2(dyjp * dxj - dxjp * dyj, dxjp * dxj + dyjp * dyj);
      statement_100:
      ctimtj = costhe(i) * costhe(j) + sinthe(i) * sinthe(j);
      stimtj = sinthe(i) * costhe(j) - sinthe(j) * costhe(i);
      aa = pi2inv * (ftan * ctimtj + flog * stimtj);
      b = pi2inv * (flog * ctimtj - ftan * stimtj);
      vtang = vtang - b * q(j) + gamma * aa;
    }
    cp(i) = 1.f - vtang * vtang;
    //C      write(6,1010)xmid,cp(i)
  }
  //C
}

//C
//C**********************************************************
//C
void
fandm(
  common& cmn,
  float const& sinalf,
  float const& cosalf,
  float& cl)
{
  common_write write(cmn);
  // COMMON bod
  arr_cref<float> x(cmn.x, dimension(100));
  arr_cref<float> y(cmn.y, dimension(100));
  // COMMON cpd
  arr_cref<float> cp(cmn.cp, dimension(100));
  //
  //C
  //C  compute and print out cd,cl,cmle
  //C
  float cfx = 0.0f;
  float cfy = 0.0f;
  float cm = 0.0f;
  //C
  int i = fem::int0;
  float xmid = fem::float0;
  float ymid = fem::float0;
  float dx = fem::float0;
  float dy = fem::float0;
  FEM_DO_SAFE(i, 1, cmn.nodtot) {
    xmid = .5f * (x(i) + x(i + 1));
    ymid = .5f * (y(i) + y(i + 1));
    dx = x(i + 1) - x(i);
    dy = y(i + 1) - y(i);
    cfx += cp(i) * dy;
    cfy = cfy - cp(i) * dx;
    cm += cp(i) * (dx * xmid + dy * ymid);
  }
  float cd = cfx * cosalf + cfy * sinalf;
  cl = cfy * cosalf - cfx * sinalf;
  /*
  write(6,
    "(/,/,/,/,'    cd =',f8.5,'    cl =',f8.5,'    cm =',f8.5)"), cd,
    cl, cm;
    */
  //C
}

struct gauss_save
{
  fem::variant_bindings cof_bindings;
};

//C
//C**********************************************************
//C
void
gauss(
  common& cmn,
  int const& nrhs)
{
  FEM_CMN_SVE(gauss);
  common_variant cof(cmn.common_cof, sve.cof_bindings);
  if (is_called_first_time) {
    using fem::mbr; // member of variant common or equivalence
    {
      mbr<float> a(dimension(101, 111));
      mbr<int> neqns;
      cof.allocate(), a, neqns;
    }
  }
  arr_ref<float, 2> a(cof.bind<float>(), dimension(101, 111));
  int const& neqns = cof.bind<int>();
  int np = fem::int0;
  int ntot = fem::int0;
  int i = fem::int0;
  int im = fem::int0;
  int imax = fem::int0;
  float amax = fem::float0;
  int j = fem::int0;
  float temp = fem::float0;
  float r = fem::float0;
  int k = fem::int0;
  int l = fem::int0;
  int ip = fem::int0;
  //C
  //C  solution of linear algebraic system by
  //C  gaussian elimination with partial pivoting
  //C
  //C          [a] = coefficient matrix
  //C          neqns = number of equations
  //C          nrhs = number of right-hand sides
  //C
  //C          right-hand sides and solutions stored in
  //C          columns neqns+1 thru neqns+nrhs of a
  //C
  np = neqns + 1;
  ntot = neqns + nrhs;
  //C
  //C  gauss reduction
  //C
  FEM_DO_SAFE(i, 2, neqns) {
    //C
    //C  search for largest entry in (i-1)th column
    //C  on or below main diagonal
    //C
    im = i - 1;
    imax = im;
    amax = fem::abs(a(im, im));
    FEM_DO_SAFE(j, i, neqns) {
      if (amax >= fem::abs(a(j, im))) {
        goto statement_110;
      }
      imax = j;
      amax = fem::abs(a(j, im));
      statement_110:;
    }
    //C
    //C  switch (i-1)th and imaxth equations
    //C
    if (imax != im) {
      goto statement_140;
    }
    FEM_DO_SAFE(j, im, ntot) {
      temp = a(im, j);
      a(im, j) = a(imax, j);
      a(imax, j) = temp;
    }
    //C
    //C  eliminate (i-1)th unknown from
    //C  ith thru neqnsth equations
    //C
    statement_140:
    FEM_DO_SAFE(j, i, neqns) {
      r = a(j, im) / a(im, im);
      FEM_DO_SAFE(k, i, ntot) {
        a(j, k) = a(j, k) - r * a(im, k);
      }
    }
  }
  //C
  //C  back substitution
  //C
  FEM_DO_SAFE(k, np, ntot) {
    a(neqns, k) = a(neqns, k) / a(neqns, neqns);
    FEM_DO_SAFE(l, 2, neqns) {
      i = neqns + 1 - l;
      ip = i + 1;
      FEM_DO_SAFE(j, ip, neqns) {
        a(i, k) = a(i, k) - a(i, j) * a(j, k);
      }
      a(i, k) = a(i, k) / a(i, i);
    }
  }
  //C
}

//C
//C**********************************************************
//C
void
indata(
  common& cmn)
{
  // COMMON par
  int& naca = cmn.naca;
  float& epsmax = cmn.epsmax;
  float& ptmax = cmn.ptmax;
  //
  //C
  //C  set parameters of body shape, flow
  //C  situation, and node distribution
  //C
  //C  user must input:
  //C  nlower = number of nodes on lower surface
  //C  nupper = number of nodes on upper surface
  //C  plus data on body
  //C
  //C      write(6,*)'input nlower,nupper'
  //C      read(5,*)nlower,nupper
  //C      write(6,*)'input naca number'
  //C      read(5,*)naca
  //C
  cmn.nlower = 50;
  cmn.nupper = 50;
  naca = 2412;
  int ieps = naca / 1000;
  int iptmax = naca / 100 - 10 * ieps;
  int itau = naca - 1000 * ieps - 100 * iptmax;
  epsmax = ieps * 0.01f;
  ptmax = iptmax * 0.1f;
  cmn.tau = itau * 0.01f;
  if (ieps < 10) {
    return;
  }
  ptmax = 0.2025f;
  epsmax = 2.6595f * fem::pow3(ptmax);
  //C
}

void
program_panel(
  int argc,
  char const* argv[])
{
  common cmn(argc, argv);
  common_write write(cmn);
  // COMMON commonymous
  arr_ref<float> t_lift(cmn.t_lift, dimension(10000));
  //
  //C
  //C        smith-hess panel method for single
  //C        element lifting airfoil in 2-d
  //C        incompressible flow
  //C
  float pi = 3.1415926585f;
  //C
  indata(cmn);
  setup(cmn);
  //C
  float alpha = 5;
  float twist = 5;
  float c_tip = 3;
  float c_root = 1;
  float s_span = 5;
  //C
  float q_dyn = 1.225f * 0.5f * fem::pow2(5);
  alpha += twist;
  float d_s = s_span / 10000.f;
  float d_twist = twist / 10000.f;
  float d_chord = (c_root - c_tip) / 10000.f;
  //C
  float tlift = 0.f;
  //C
  int i = fem::int0;
  float cosalf = fem::float0;
  float sinalf = fem::float0;
  float cl = fem::float0;
  float chord = fem::float0;
  float area = fem::float0;
  FEM_DOSTEP(i, 1, 10000, 1) {
    //C
    cosalf = fem::cos(alpha * pi / 180.f);
    sinalf = fem::sin(alpha * pi / 180.f);
    cofish(cmn, sinalf, cosalf);
    gauss(cmn, 1);
    veldis(cmn, sinalf, cosalf);
    fandm(cmn, sinalf, cosalf, cl);
    //C
    alpha = alpha - d_twist;
    chord = c_root - d_chord * i;
    area = d_s * chord;
    //C
    t_lift(i) = cl * q_dyn * area;
    //C
  }
  //C
  FEM_DOSTEP(i, 1, 10000, 1) {
    tlift += t_lift(i);
  }
  //C
  //C        print*, "Lift coef:",2.*t_lift/(q_dyn*(c_root+c_tip)*0.5*s_span)
  write(6, star), "Total Lift :", 2.f * tlift, "  [N]";
  //C      go to 100
  //C
  FEM_STOP(0);
}

} // namespace placeholder_please_replace

int
main(
  int argc,
  char const* argv[])
{
  return fem::main_with_catch(
    argc, argv,
    hw::program_panel);
}
