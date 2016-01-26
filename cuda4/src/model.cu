#include <fem.hpp> // Fortran EMulation library of fable module

#include "constants.h"
#include "arrays.h"

#include <ppl.h>
#include <chrono>

#define FEM_DO_SAFE_1(x,a,b) for ( int32_t x = a; x<=b; ++x )
#define FEM_DOSTEP_1(x, a, b, c) for (int32_t x = a; x <= b; x = x + c)

namespace hw {

typedef array_1d_mn_cpu<float, constants::n > data_arr_float;

data_arr_float costhe;
data_arr_float sinthe;


struct context
{
    float cof_memory[101][111];
    array_2d_mn<float, 101, 111> cof;

    data_arr_float cp;

    context() :
    cof(&cof_memory[0][0])
    {

    }
};


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

common_par  par;

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



struct common_commonymous
{
  arr<float> t_lift;

  common_commonymous() :
    t_lift(dimension(10000), fem::fill0)
  {}
};

struct common :
  fem::common,
  common_bod,
  common_commonymous
{
  fem::cmn_sve cofish_sve;
  fem::cmn_sve veldis_sve;
  fem::cmn_sve gauss_sve;

  float cof_memory[101][111];
  array_2d_mn<float, 101, 111> cof;

  common(
    int argc,
    char const* argv[])
  :
    fem::common(argc, argv)
    ,cof(&cof_memory[0][0])
  {}
};


inline int kutta(const common& c)
{
    return c.nodtot + 1;
}

inline int neqns(const common& c)
{
    return kutta(c);
}

//C
//C**********************************************************
//C
void
naca45(
  const common& cmn,
  const common_par& par,
  float const& z,
  float& thick,
  float& camber,
  float& beta)
{
  float epsmax = par.epsmax;
  float ptmax = par.ptmax;
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
  thick = 5.f * par.tau * (.2969f * fem::sqrt(z) - z * (.126f + z * (
    .3537f - z * (.2843f - z * .1015f))));
  statement_100:
  if (epsmax == 0.f) {
    goto statement_130;
  }
  if (par.naca > 9999) {
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
  const common& cmn,
  const common_par& par,
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
  naca45(cmn, par, z, thick, camber, beta);
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
  FEM_DO_SAFE_1(nsurf, 1, 2) {
    FEM_DO_SAFE_1(n, 1, npoints) {
      fract = fem::ffloat(n - 1) / fem::ffloat(npoints);
      z = .5f * (1.f - fem::cos(pi * fract));
      i = nstart + n;
      body(cmn, par, z, sign, x(i), y(i));
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
  FEM_DO_SAFE_1(i, 1, nodtot) {
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
  const common& cmn,
  float sinalf,
  float cosalf,
  context& c    )
{
  //FEM_CMN_SVE(cofish);
  int nodtot = cmn.nodtot;
  arr_cref<float> x(cmn.x, dimension(100));
  arr_cref<float> y(cmn.y, dimension(100));

  const float  pi2inv = constant_functions::pi2inv();

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
  auto kutta = hw::kutta(cmn);
  //C
  //C  initialize coefficients
  //C
  FEM_DO_SAFE_1(j, 1, kutta) {
   c.cof(kutta, j) = 0.0f;
  }
  //C
  //C  set vn=0. at midpoint of ith panel
  //C
  FEM_DO_SAFE_1(i, 1, nodtot) {
    xmid = .5f * (x(i) + x(i + 1));
    ymid = .5f * (y(i) + y(i + 1));
    c.cof(i, kutta) = 0.0f;
    //C
    //C  find contribution of jth panel
    //C
    FEM_DO_SAFE_1(j, 1, nodtot) {
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
      c.cof(i, j) = pi2inv * (ftan * ctimtj + flog * stimtj);
      b = pi2inv * (flog * ctimtj - ftan * stimtj);
      c.cof(i, kutta) += b;
      if ((i > 1) && (i < nodtot)) {
        goto statement_110;
      }
      //C
      //C  if ith panel touches trailing edge, add contribution
      //C    to kutta condition
      //C
      c.cof(kutta, j) = c.cof(kutta, j) - b;
      c.cof(kutta, kutta) += c.cof(i, j);
      statement_110:;
    }
    //C
    //C  fill in known sides
    //C
    c.cof(i, kutta + 1) = sinthe(i) * cosalf - costhe(i) * sinalf;
  }
  c.cof(kutta, kutta + 1) = -(costhe(1) + costhe(nodtot)) * cosalf - (
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
  context& c,
  float const& sinalf,
  float const& cosalf)
{
  //FEM_CMN_SVE(veldis);
  const int& nodtot = cmn.nodtot;
  arr_cref<float> x(cmn.x, dimension(100));
  arr_cref<float> y(cmn.y, dimension(100));

  //arr_ref<float> cp(cmn.cp, dimension(100));
  const float pi2inv = constant_functions::pi2inv();

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

  auto kutta = hw::kutta(cmn);
  //C
  //C  compute and print out pressure distribution
  //C
  //C      write(6,1000)
  //C
  //C  retrieve solution from a-matrix
  //C
  FEM_DO_SAFE_1(i, 1, nodtot) {
    q(i) = c.cof(i, kutta + 1);
  }
  gamma = c.cof(kutta, kutta + 1);
  //C
  //C  find vt and cp at midpoint of ith panel
  //C
  FEM_DO_SAFE_1(i, 1, nodtot) {
    xmid = .5f * (x(i) + x(i + 1));
    ymid = .5f * (y(i) + y(i + 1));
    vtang = cosalf * costhe(i) + sinalf * sinthe(i);
    //C
    //C  add contributions of jth panel
    //C
    FEM_DO_SAFE_1(j, 1, nodtot) {
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
    c.cp(i) = 1.f - vtang * vtang;
    //C      write(6,1010)xmid,cp(i)
  }
  //C
}

//C
//C**********************************************************
//C
void
fandm(
  const common& cmn,
  const context& c,
  float const& sinalf,
  float const& cosalf,
  float& cl)
{
  arr_cref<float> x(cmn.x, dimension(100));
  arr_cref<float> y(cmn.y, dimension(100));
  // COMMON cpd
//  arr_cref<float> cp(cmn.cp, dimension(100));
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
  FEM_DO_SAFE_1(i, 1, cmn.nodtot) {
    xmid = .5f * (x(i) + x(i + 1));
    ymid = .5f * (y(i) + y(i + 1));
    dx = x(i + 1) - x(i);
    dy = y(i + 1) - y(i);
    cfx += c.cp(i) * dy;
    cfy = cfy - c.cp(i) * dx;
    cm += c.cp(i) * (dx * xmid + dy * ymid);
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
    const common& cmn,
    context& c,
    int nrhs)
{
  //FEM_CMN_SVE(gauss);
  
 
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
  auto neqns = hw::neqns(cmn);
  np = neqns + 1;
  ntot = neqns + nrhs;
  //C
  //C  gauss reduction
  //C
  FEM_DO_SAFE_1(i, 2, neqns) {
    //C
    //C  search for largest entry in (i-1)th column
    //C  on or below main diagonal
    //C
    im = i - 1;
    imax = im;
    amax = fem::abs(c.cof(im, im));
    FEM_DO_SAFE_1(j, i, neqns) {
      if (amax >= fem::abs(c.cof(j, im))) {
        goto statement_110;
      }
      imax = j;
      amax = fem::abs(c.cof(j, im));
      statement_110:;
    }
    //C
    //C  switch (i-1)th and imaxth equations
    //C
    if (imax != im) {
      goto statement_140;
    }
    FEM_DO_SAFE_1(j, im, ntot) {
      temp = c.cof(im, j);
      c.cof(im, j) = c.cof(imax, j);
      c.cof(imax, j) = temp;
    }
    //C
    //C  eliminate (i-1)th unknown from
    //C  ith thru neqnsth equations
    //C
    statement_140:
    FEM_DO_SAFE_1(j, i, neqns) {
      r = c.cof(j, im) / c.cof(im, im);
      FEM_DO_SAFE_1(k, i, ntot) {
          c.cof(j, k) = c.cof(j, k) - r * c.cof(im, k);
      }
    }
  }
  //C
  //C  back substitution
  //C
  FEM_DO_SAFE_1(k, np, ntot) {
      c.cof(neqns, k) = c.cof(neqns, k) / c.cof(neqns, neqns);
    FEM_DO_SAFE_1(l, 2, neqns) {
      i = neqns + 1 - l;
      ip = i + 1;
      FEM_DO_SAFE_1(j, ip, neqns) {
          c.cof(i, k) = c.cof(i, k) - c.cof(i, j) * c.cof(j, k);
      }
      c.cof(i, k) = c.cof(i, k) / c.cof(i, i);
    }
  }
  //C
}

//C
//C**********************************************************
//C
void
indata(
  common& cmn,
  common_par& par
    )
{
  // COMMON par
  int& naca = par.naca;
  float& epsmax = par.epsmax;
  float& ptmax = par.ptmax;
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
  par.tau = itau * 0.01f;
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
  indata(cmn, par);
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

  std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
  concurrency::parallel_for(1, 10000, [ &t_lift, c_root, d_chord, d_s, alpha, pi, d_twist, q_dyn, &cmn ] ( int i )
  {
          context c;

          float chord = fem::float0;
          float area = fem::float0;
          float cl = fem::float0;
          //C
          auto cosalf = fem::cos( (alpha - (d_twist * (i - 1)))  * pi / 180.f);
          auto sinalf = fem::sin( (alpha - (d_twist * (i - 1)))  * pi / 180.f);
          cofish(cmn, sinalf, cosalf, c);
          gauss(cmn, c, 1);
          veldis(cmn, c, sinalf, cosalf);
          fandm(cmn, c, sinalf, cosalf, cl);
          chord = c_root - d_chord * i;
          area = d_s * chord;
          t_lift(i) = cl * q_dyn * area;
  });

  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

  std::cout << "Filtering on host took: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms" << std::endl;

  //C
  FEM_DOSTEP_1(i, 1, 10000, 1) {
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