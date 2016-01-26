#include <fem.hpp> // Fortran EMulation library of fable module

#include "constants.h"
#include "arrays.h"


#include <ppl.h>
#include <chrono>

#include <thrust/for_each.h>
#include "cuda_helper.h"

#define FEM_DO_SAFE_1(x,a,b) for ( int32_t x = a; x<=b; ++x )
#define FEM_DOSTEP_1(x, a, b, c) for (int32_t x = a; x <= b; x = x + c)

namespace hw {

struct context
{
    array_2d_mn_fixed<float, 101, 111> cof;
    array_1d_mn_fixed<float, 100> cp;
};


using namespace fem::major_types;

__host__ __device__ inline int zero_int()
{
    return 0;
}

__host__ __device__ inline float zero_float()
{
    return 0.0f;
}


struct common_par
{
  int naca;
  float tau;
  float epsmax;
  float ptmax;

  common_par() :
    naca(zero_int()),
    tau(zero_float()),
    epsmax(zero_float()),
    ptmax(zero_float())
  {}
};



struct common_input
{
  int nlower;
  int nupper;
  int nodtot;

  array_1d_mn_fixed<float, 101> x;
  array_1d_mn_fixed<float, 101> y;

  array_1d_mn_fixed<float, 101> costhe;
  array_1d_mn_fixed<float, 101> sinthe;

  common_input() :
    nlower(zero_int()),
    nupper(zero_int()),
    nodtot(zero_int())
  {

  }

  void patch_pointers()
  {

  }

};

struct common_output
{
  array_1d_mn_fixed<float, 10001> t_lift;

  common_output() 
  {}
};

struct common :
  common_input,
  common_output
{
  array_2d_mn_fixed<float, 101, 111> cof;

  common(
    int argc,
    char const* argv[])
  {}
};


__host__ __device__ inline int kutta(const common& c)
{
    return c.nodtot + 1;
}

__host__ __device__ inline int neqns(const common& c)
{
    return kutta(c);
}

//C
//C**********************************************************
//C
__host__ __device__ void
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
  float dcamdx = zero_float();
  float w = zero_float();
  //C
  //C  evaluate thickness and camber
  //C  for naca 4- or 5-digit airfoil
  //C
  thick = 0.0f;
  if (z < 1.e-10f) {
    goto statement_100;
  }
  thick = 5.f * par.tau * (.2969f * sqrt(z) - z * (.126f + z * (
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
  beta = atan(dcamdx);
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
  float thick = zero_float();
  float camber = zero_float();
  float beta = zero_float();
  naca45(cmn, par, z, thick, camber, beta);
  x = z - sign * thick * sin(beta);
  y = camber + sign * thick * cos(beta);
  //C
}

//C
//C**********************************************************
//C
void
setup( const common_par& par, common& cmn    )
{
  fem::common  w;
  common_write write(w);
  // COMMON bod
  int& nlower = cmn.nlower;
  int& nupper = cmn.nupper;
  int& nodtot = cmn.nodtot;

  auto& x = cmn.x;
  auto& y = cmn.y;
  auto& costhe = cmn.costhe;
  auto& sinthe = cmn.sinthe;

  const float pi = constant_functions::pi();
  
  //C
  //C  set coordinates of nodes on body surface
  //C
  write(6, "(/,/,/,' body shape',/,/,4x,'x',9x,'y',/)");
  int npoints = nlower;
  float sign = -1.0f;
  int nstart = 0;
  int n = zero_int();
  float fract = zero_float();
  float z = zero_float();
  int i = zero_int();
  FEM_DO_SAFE_1(nsurf, 1, 2) {
    FEM_DO_SAFE_1(n, 1, npoints) {
      fract = fem::ffloat(n - 1) / fem::ffloat(npoints);
      z = .5f * (1.f - cos(pi * fract));
      i = nstart + n;
      float xa;
      float ya;

      body(cmn, par, z, sign, xa, ya);
      x(i) = xa;
      y(i) = ya;

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
  float dx = zero_float();
  float dy = zero_float();
  float dist = zero_float();
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
__host__ __device__ void
cofish(
  const common& cmn,
  float sinalf,
  float cosalf,
  context& c    )
{
  //FEM_CMN_SVE(cofish);
  int nodtot = cmn.nodtot;
  const auto& x = cmn.x;
  const auto& y = cmn.y;
  const auto& costhe = cmn.costhe;
  const auto& sinthe = cmn.sinthe;

  const float  pi2inv = constant_functions::pi2inv();

  float xmid = zero_float();
  float ymid = zero_float();
  float flog = zero_float();
  float ftan = zero_float();
  float dxj = zero_float();
  float dxjp = zero_float();
  float dyj = zero_float();
  float dyjp = zero_float();
  float ctimtj = zero_float();
  float stimtj = zero_float();
  float b = zero_float();
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
      flog = .5f * log((dxjp * dxjp + dyjp * dyjp) / (dxj *
        dxj + dyj * dyj));
      ftan = atan2(dyjp * dxj - dxjp * dyj, dxjp * dxj + dyjp * dyj);
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
__host__ __device__ void
veldis(
  const common& cmn,
  context& c,
  float const& sinalf,
  float const& cosalf)
{
  //FEM_CMN_SVE(veldis);
  const int& nodtot = cmn.nodtot;

  auto& x = cmn.x;
  auto& y = cmn.y;

  const auto& costhe = cmn.costhe;
  const auto& sinthe = cmn.sinthe;

  //arr_ref<float> cp(cmn.cp, dimension(100));
  const float pi2inv = constant_functions::pi2inv();

  int i = zero_int();
  array_1d_mn_fixed<float, 150> q;
  float gamma = zero_float();
  float xmid = zero_float();
  float ymid = zero_float();
  float vtang = zero_float();
  int j = zero_int();
  float flog = zero_float();
  float ftan = zero_float();
  float dxj = zero_float();
  float dxjp = zero_float();
  float dyj = zero_float();
  float dyjp = zero_float();
  float ctimtj = zero_float();
  float stimtj = zero_float();
  float aa = zero_float();
  float b = zero_float();

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
      flog = .5f * log((dxjp * dxjp + dyjp * dyjp) / (dxj *
        dxj + dyj * dyj));
      ftan = atan2(dyjp * dxj - dxjp * dyj, dxjp * dxj + dyjp * dyj);
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
__host__ __device__ void
fandm(
  const common& cmn,
  const context& c,
  float const& sinalf,
  float const& cosalf,
  float& cl)
{
    const auto& x = cmn.x;
    const auto& y = cmn.y;
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
  int i = zero_int();
  float xmid = zero_float();
  float ymid = zero_float();
  float dx = zero_float();
  float dy = zero_float();
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
__host__ __device__ void
gauss(
    const common& cmn,
    context& c,
    int nrhs)
{
  int np = zero_int();
  int ntot = zero_int();
  int i = zero_int();
  int im = zero_int();
  int imax = zero_int();
  float amax = zero_float();
  float temp = zero_float();
  float r = zero_float();
  int ip = zero_int();
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
      if (amax >= abs(c.cof(j, im))) {
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

struct kernel
{
    common* m_c;
    float   m_c_root;
    float   m_d_chord;
    float   m_d_s;
    float   m_alpha;
    float   m_pi;
    float   m_d_twist;
    float   m_q_dyn;

    __host__ kernel(common* c, float c_root, float d_chord, float d_s, float alpha, float pi, float d_twist, float q_dyn ) : m_c(c)
    { 
        m_c_root = c_root;
        m_d_chord = d_chord;
        m_d_s = d_s;
        m_alpha = alpha;
        m_pi = pi;
        m_d_twist = d_twist;
        m_q_dyn = q_dyn;
    }

    __device__ void operator() (int32_t i)
    {
        i = i + 1;
        common& cmn = *m_c;

        float c_root = m_c_root;
        float d_chord = m_d_chord;
        float d_s = m_d_s;
        float alpha = m_alpha;
        float pi = m_pi;
        float d_twist = m_d_twist;
        float q_dyn = m_q_dyn;

        context c;

        float chord = zero_float();
        float area = zero_float();
        float cl = zero_float();
        //C
        auto cosalf = cos((alpha - (d_twist * (i - 1)))  * pi / 180.f);
        auto sinalf = sin((alpha - (d_twist * (i - 1)))  * pi / 180.f);
        cofish(cmn, sinalf, cosalf, c);
        gauss(cmn, c, 1);
        veldis(cmn, c, sinalf, cosalf);
        fandm(cmn, c, sinalf, cosalf, cl);
        chord = c_root - d_chord * i;
        area = d_s * chord;
        cmn.t_lift(i) = cl * q_dyn * area;
    }
};


static inline common* allocate_device_common_buffer(const common& c)
{
    common* pointer;
    size_t  s = sizeof(common);

    ::cuda::throw_if_failed(cudaMalloc(&pointer, s));
    ::cuda::throw_if_failed(cudaMemset(pointer, 0, s));
    ::cuda::throw_if_failed(cudaMemcpy(pointer, &c, s, cudaMemcpyHostToDevice));

    return pointer;
}


void launch_kernel(const common& c, float c_root, float d_chord, float d_s, float alpha, float pi, float d_twist, float q_dyn )
{
    auto pointer = allocate_device_common_buffer(c);

    kernel k(pointer, c_root, d_chord, d_s, alpha, pi, d_twist, q_dyn);

    auto cb = thrust::make_counting_iterator<std::uint32_t>(0);
    auto ce = cb + 10000 ;

    thrust::for_each(cb, ce, kernel(pointer, c_root, d_chord, d_s, alpha, pi, d_twist, q_dyn) );

    ::cuda::throw_if_failed(cudaDeviceSynchronize());
}


void
program_panel(
  int argc,
  char const* argv[])
{
  common cmn(argc, argv);
  common_par par;
  fem::common  w;
  common_write write(w);

  
  auto& t_lift = cmn.t_lift;
  //
  //C
  //C        smith-hess panel method for single
  //C        element lifting airfoil in 2-d
  //C        incompressible flow
  //C
  float pi = 3.1415926585f;
  //C
  indata(cmn, par);
  setup(par, cmn);
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
  
  //launch_kernel(cmn, c_root, d_chord, d_s, alpha, pi, d_twist, q_dyn);

  std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
  concurrency::parallel_for(1, 10000, [ &t_lift, c_root, d_chord, d_s, alpha, pi, d_twist, q_dyn, &cmn ] ( int i )
  {
          context c;

          float chord = zero_float();
          float area = zero_float();
          float cl = zero_float();
          //C
          auto cosalf = cos( (alpha - (d_twist * (i - 1)))  * pi / 180.f);
          auto sinalf = sin( (alpha - (d_twist * (i - 1)))  * pi / 180.f);
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
