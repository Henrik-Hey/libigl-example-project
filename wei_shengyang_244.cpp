#include <igl/PI.h>
#include <igl/upsample.h>

// Wei Shengyang Lifted Loop WT: 2.4.4

// Common/Helper Functions 

double WT_scale_common(
    const int n
) {
    double p_1 = std::cos((2. * M_PI) / n);
    double p_2 = 0.375 + (0.25*p_1); 
    return p_2*p_2;
}

/**
 * Lifting I. 
 * The first lifting step subtracts a filtered version of the boundary vertices in Vnew from those in Vold.
 * Let v be a boundary vertex in Vold, and v1 and v2 be boundary vertices from Vnew. The relative positions of v1
 * and v2 to v are specified by the mask shown in Figure 2.16(a). Then, the updated value v` of v is given by
*/

void WT_Lifting_1(
    const Eigen::Vector3d v,
    const Eigen::Vector3d v1,
    const Eigen::Vector3d v2,
    Eigen::Vector3d& v_prime
) {
    v_prime = v - (0.25 * (v1 + v2));
}

void WT_Lifting_1_inverse(
    const Eigen::Vector3d v_prime,
    const Eigen::Vector3d v1,
    const Eigen::Vector3d v2,
    Eigen::Vector3d& v
) {
    v = v_prime + (0.25 * (v1 + v2));
}

/**
 * Lifting II. 
 * The second lifting step subtracts a filtered version of the boundary vertices in Vold from those in
 * Vnew. Let v be a boundary vertex in Vnew, and v1 and v2 be vertices from Vold. The relative positions of v1 and v2
 * to v are specified by the mask shown in Figure 2.16(b). Then, the updated value v` of v is given by
*/

void WT_Lifting_2(
    const Eigen::Vector3d v,
    const Eigen::Vector3d v1,
    const Eigen::Vector3d v2,
    Eigen::Vector3d& v_prime
) {
    v_prime = v - (0.5 * (v1 + v2));
}

void WT_Lifting_2_inverse(
    const Eigen::Vector3d v_prime,
    const Eigen::Vector3d v1,
    const Eigen::Vector3d v2,
    Eigen::Vector3d& v
) {
    v = v_prime + (0.5 * (v1 + v2));
}

// Scalar delta function for Lifting III (3)

double scalar_delta(
    const int n
) {
    double n_1 = 1 / n;
    double p_4 = 1 - ((8./5.) * WT_scale_common(n));
    return n_1 * p_4;
}

/**
 * Lifting III. 
 * The third lifting step subtracts a filtered version of the interior vertices in Vnew from those in Vold.
 * Let v be a vertex in Vold, and v1, v2 ..., and vn be the n 1-ring neighbours of v, which are also in Vnew. The
 * relative positions of v1, v2 ..., and vn to v are specified by the mask shown in Figure 2.16(c). Then, the updated
 * value v` of v is given by
*/

Eigen::Vector3d WT_Lifting_3_getScalarDeltaSum(
    const std::vector<Eigen::Vector3d> vertices
) {
    Eigen::Vector3d sum;
    sum << 0, 0, 0;

    for(
        std::vector<Eigen::Vector3d>::const_iterator it = vertices.begin();
        it != vertices.end();
        it++
    ){
        sum += *it;
    }

    return scalar_delta(vertices.size()) * sum;
}

void WT_Lifting_3(
    const std::vector<Eigen::Vector3d> vertices,
    const Eigen::Vector3d v,
    Eigen::Vector3d& v_prime
) {
    v_prime = v - WT_Lifting_3_getScalarDeltaSum(vertices);
}

void WT_Lifting_3_inverse(
    const Eigen::Vector3d v_prime,
    const std::vector<Eigen::Vector3d> vertices,
    Eigen::Vector3d& v
) {
    v = v_prime + WT_Lifting_3_getScalarDeltaSum(vertices);
}

/**
 * Scaling. 
 * The scaling step divides the interior vertices in Vold by a scalar. Let v be an interior vertex of valence
 * n in Vold. The updated value of v` of v is given by
*/

void WT_Scaling_inverse(
    const Eigen::Vector3d v_prime, 
    const double n,
    Eigen::Vector3d& v
) {
    v = v_prime * ((8 / 5) * WT_scale_common(n));    
}

/**
 * Lifting IV. 
 * The fourth lifting step subtracts a filtered version of the vertices in Vold from the interior vertices in
 * Vnew. Let v be an interior vertex in Vnew, and v1, v2, v3, and v4 be vertices from Vold. The relative positions of v1,
*/

void WT_Lifting_4(
    const Eigen::Vector3d v,
    const Eigen::Vector3d v1,
    const Eigen::Vector3d v2,
    const Eigen::Vector3d v3,
    const Eigen::Vector3d v4,
    Eigen::Vector3d& v_prime
) {
    const Eigen::Vector3d three_eigths = (3./8.) * (v1 + v2);
    const Eigen::Vector3d one_eigth = (1./8.) * (v3 + v4); 

    v_prime = v - (three_eigths + one_eigth); 
}

void WT_Lifting_4_inverse(
    Eigen::Vector3d& v_prime,
    const Eigen::Vector3d v1,
    const Eigen::Vector3d v2,
    const Eigen::Vector3d v3,
    const Eigen::Vector3d v4,
    Eigen::Vector3d& v
) {
    const Eigen::Vector3d three_eigths = (3./8.) * (v1 + v2);
    const Eigen::Vector3d one_eigth = (1./8.) * (v3 + v4); 

    v = v_prime + (three_eigths + one_eigth); 
}


/** 
 * Lifting V. 
 * The fifth lifting step subtracts a filtered version of the boundary vertices in Vnew from those in Vold
 * Let v1, v2, v3, and v4 be boundary vertices in Vold, and v be a boundary vertex from Vnew. The relative positions
 * of v1, v2, v3, and v4 to v are specified by the mask shown in Figure 2.16(e). Then, the updated values of v`1, v`2, v`3,
 * and v`4 of v1, v2, v3, and v4 are given by
*/

void WT_Lifting_5(
    const Eigen::Vector3d v,
    const Eigen::Vector3d v1,
    const Eigen::Vector3d v2,
    const Eigen::Vector3d v3,
    const Eigen::Vector3d v4,
    Eigen::Vector3d& v1_prime,
    Eigen::Vector3d& v2_prime,
    Eigen::Vector3d& v3_prime,
    Eigen::Vector3d& v4_prime
) {
    const double n_14 = -0.525336;
    const double n_23 = 0.189068;

    v1_prime = v1 - (n_14 * v);
    v2_prime = v2 - (n_23 * v);
    v3_prime = v3 - (n_23 * v);
    v4_prime = v4 - (n_14 * v);
}

void WT_Lifting_5_inverse(
    const Eigen::Vector3d v,
    const Eigen::Vector3d v1_prime,
    const Eigen::Vector3d v2_prime,
    const Eigen::Vector3d v3_prime,
    const Eigen::Vector3d v4_prime,
    Eigen::Vector3d& v1,
    Eigen::Vector3d& v2,
    Eigen::Vector3d& v3,
    Eigen::Vector3d& v4
) {
    const double n_14 = -0.525336;
    const double n_23 = 0.189068;

    v1 = v1_prime + (n_14 * v);
    v2 = v2_prime + (n_23 * v);
    v3 = v3_prime + (n_23 * v);
    v4 = v4_prime + (n_14 * v);
}

/**
 * Lifting VI
*/

double WT_Coefficient_Alpha(
    const double n
) {
  return (0.375) + WT_scale_common(n); 
}

void WT_Scaling(
    const Eigen::Vector3d v,
    const double n,
    Eigen::Vector3d& v_prime 
) {
    // if(n==0) return;
    // double scale = WT_scale_common(n);
    // assert(scale != 0);

    double beta = (WT_Coefficient_Alpha(n) - 0.375) / 0.625;
    // v_prime = v / ((8. / 5.) * scale);    
    v_prime = v / beta;    
}

double WT_Coefficient_Beta(
    const double n
) {
    return (1.6) * WT_scale_common(n); 
}

double WT_Coefficient_Gamma(
    const double n
) {
    return (1. / n) * ((5./8.) - WT_scale_common(n)); 
}

double WT_Coefficient_Delta(
    const double n
) {
    return (1./n) * (1. - WT_Coefficient_Beta(n)); 
}

void WT_Get_A(
    const double n_0,
    const double n_1,
    const double n_2,
    const double n_3,
    Eigen::Matrix4d& A
) {
    // Coefficients

    const double alpha_0 = WT_Coefficient_Alpha(n_0);
    const double alpha_1 = WT_Coefficient_Alpha(n_1);
    const double alpha_2 = WT_Coefficient_Alpha(n_2);
    const double alpha_3 = WT_Coefficient_Alpha(n_3);

    const double gamma_0 = WT_Coefficient_Gamma(n_0);
    const double gamma_1 = WT_Coefficient_Gamma(n_1);
    const double gamma_2 = WT_Coefficient_Gamma(n_2);
    const double gamma_3 = WT_Coefficient_Gamma(n_3);

    /** START OF VARIABLE INITIALIZATION **/

    /** 
     * IDENTITY BOYOS
    */
    double a_00 = 0;
    double a_11 = 0;
    double a_22 = 0;
    double a_33 = 0;

    /** 
     * Mirror-able
     * These values can be mirrored. i.e. a_01 == a_21, a_02 == a_31, etc...
    */
    double a_01 = 0;
    double a_02 = 0;
    double a_03 = 0;
    double a_12 = 0;
    double a_13 = 0;
    double a_23 = 0;

    /** END OF VARIABLE INITIALIZATION **/

    // Hair-loss-math-time (its not hard just exhaustive)
    // first start with the IDENTITY BOYOS

    a_00 += std::pow(alpha_0, 2);
    a_00 += std::pow(gamma_1, 2);
    a_00 += std::pow(gamma_2, 2);
    a_00 += std::pow(gamma_3, 2);
    a_00 += (1 / 256) * (n_0 - 3);
    a_00 += (5 / 32) * n_0;

    a_11 += std::pow(gamma_0, 2);
    a_11 += std::pow(alpha_1, 2);
    a_11 += std::pow(gamma_2, 2);
    a_11 += std::pow(gamma_3, 2);
    a_11 += (1 / 256) * (n_1 - 3);
    a_11 += (5 / 32) * n_1;

    a_22 += std::pow(gamma_0, 2);
    a_22 += std::pow(gamma_1, 2);
    a_22 += std::pow(alpha_2, 2);
    a_22 += (1 / 256) * (n_2 - 3);
    a_22 += (5 / 32) * n_2;

    a_33 += std::pow(gamma_0, 2);
    a_33 += std::pow(gamma_1, 2);
    a_33 += std::pow(alpha_3, 2);
    a_33 += (1 / 256) * (n_3 - 3);
    a_33 += (5 / 32) * n_3;

    // lastly the mirror bois

    a_01 += alpha_0 * gamma_0;
    a_01 += gamma_1 * alpha_1;
    a_01 += std::pow(gamma_2, 2);
    a_01 += std::pow(gamma_3, 2);
    a_01 += (21 / 64);

    a_02 += alpha_0 * gamma_0;
    a_02 += std::pow(gamma_1, 2);
    a_02 += gamma_2 * alpha_2;
    a_02 += (85 / 256);

    a_03 += alpha_0 * gamma_0;
    a_03 += std::pow(gamma_1, 2);
    a_03 += gamma_3 * alpha_3;
    a_03 += (85 / 256);
    
    a_12 += std::pow(gamma_0, 2);
    a_12 += alpha_1 * gamma_1;
    a_12 += gamma_2 * alpha_2;
    a_12 += (85 / 256);

    a_13 += std::pow(gamma_0, 2);
    a_13 += alpha_1 * gamma_1;
    a_13 += gamma_3 * alpha_3;
    a_13 += (85 / 256);

    a_23 += std::pow(gamma_0, 2);
    a_23 += std::pow(gamma_1, 2);
    a_23 += (1 / 64);

    A << 
        a_00, a_01, a_02, a_03,
        a_01, a_11, a_12, a_13,
        a_02, a_12, a_22, a_23,
        a_03, a_13, a_23, a_33;

    double n0d = static_cast<double>(n_0);
    double n1d = static_cast<double>(n_1);
    double n2d = static_cast<double>(n_2);
    double n3d = static_cast<double>(n_3);

    double alpha0 = WT_Coefficient_Alpha(n0d);
    double alpha1 = WT_Coefficient_Alpha(n1d);
    double alpha2 = WT_Coefficient_Alpha(n2d);
    double alpha3 = WT_Coefficient_Alpha(n3d);

    double gamma0 = WT_Coefficient_Gamma(n0d);
    double gamma1 = WT_Coefficient_Gamma(n1d);
    double gamma2 = WT_Coefficient_Gamma(n2d);
    double gamma3 = WT_Coefficient_Gamma(n3d);

    double delta0 = WT_Coefficient_Delta(n0d);
    double delta1 = WT_Coefficient_Delta(n1d);

    A(0, 0) = alpha0 * alpha0 + gamma1 * gamma1 + gamma2 * gamma2 + gamma3 * gamma3 + (n0d - 3.0) * 0.00390625 + n0d * (0.15625);
    A(1, 1) = gamma0 * gamma0 + alpha1 * alpha1 + gamma2 * gamma2 + gamma3 * gamma3 + (n1d - 3.0) * 0.00390625 + n1d * (0.15625);
    A(2, 2) = gamma0 * gamma0 + gamma1 * gamma1 + alpha2 * alpha2 + (n2d - 2.0) * 0.00390625 + n2d * (0.15625);
    A(3, 3) = gamma0 * gamma0 + gamma1 * gamma1 + alpha3 * alpha3 + (n3d - 2.0) * 0.00390625 + n3d * (0.15625);
    A(0, 1) = A(1, 0) = alpha0 * gamma0 + gamma1 * alpha1 + gamma2 * gamma2 + gamma3 * gamma3 + 0.328125;
    A(0, 2) = A(2, 0) = alpha0 * gamma0 + gamma1 * gamma1 + gamma2 * alpha2 + 0.33203125;
    A(0, 3) = A(3, 0) = alpha0 * gamma0 + gamma1 * gamma1 + gamma3 * alpha3 + 0.33203125;
    A(1, 2) = A(2, 1) = gamma0 * gamma0 + alpha1 * gamma1 + alpha2 * gamma2 + 0.33203125;
    A(1, 3) = A(3, 1) = gamma0 * gamma0 + alpha1 * gamma1 + alpha3 * gamma3 + 0.33203125;
    A(2, 3) = A(3, 2) = gamma0 * gamma0 + gamma1 * gamma1 + 0.015625;


    std::cout << "new A: " << std::endl;
    std::cout << A << std::endl;

}

void WT_GET_B(
    const double n_0,
    const double n_1,
    Eigen::Vector4d& B
) {
    const double alpha_0 = WT_Coefficient_Alpha(n_0);
    const double alpha_1 = WT_Coefficient_Alpha(n_1);

    const double gamma_0 = WT_Coefficient_Gamma(n_0);
    const double gamma_1 = WT_Coefficient_Gamma(n_1);
    
    const double delta_0 = WT_Coefficient_Delta(n_0);
    const double delta_1 = WT_Coefficient_Delta(n_1);

    double b_0 = 0;
    double b_1 = 0;
    double b_2 = 0;

    b_0 += alpha_0 * delta_0;
    b_0 += gamma_1 * delta_1;
    b_0 += (3 / 8);

    b_1 += gamma_0 * delta_0;
    b_1 += alpha_1 * delta_1;
    b_1 += (3 / 8);

    b_2 += gamma_0 * delta_0;
    b_2 += gamma_1 * delta_1;
    b_2 += (1 / 8);

    B << 
        b_0,
        b_1,
        b_2,
        b_2;

    std::cout << "Got B: " << std::endl;
    std::cout << B << std::endl;
    std::cout << "-------" << std::endl;

    double n0d = static_cast<double>(n_0);
    double n1d = static_cast<double>(n_1);

    double alpha0 = WT_Coefficient_Alpha(n0d);
    double alpha1 = WT_Coefficient_Alpha(n1d);

    double gamma0 = WT_Coefficient_Gamma(n0d);
    double gamma1 = WT_Coefficient_Gamma(n1d);

    double delta0 = WT_Coefficient_Delta(n0d);
    double delta1 = WT_Coefficient_Delta(n1d);

    B(0) = -(alpha0 * delta0 + gamma1 * delta1 + 0.375);
    B(1) = -(gamma0 * delta0 + alpha1 * delta1 + 0.375);
    B(2) = B(3) = -(gamma0 * delta0 + gamma1 * delta1 + 0.125);

    std::cout << "Got New B: " << std::endl;
    std::cout << B << std::endl;
    std::cout << "-------" << std::endl;
}

void WT_Solve_Weights(
    const double n_0,
    const double n_1,
    const double n_2,
    const double n_3,
    Eigen::Vector4d& W
) {
    Eigen::Matrix4d A;
    WT_Get_A(n_0, n_1, n_2, n_3, A);

    Eigen::Vector4d B;
    WT_GET_B(n_0, n_1, B);

    std::cout << "Computing A and B" << std::endl;
    std::cout << "valence v0 :" << n_0 << std::endl;
    std::cout << "valence v1 :" << n_1 << std::endl;
    std::cout << "valence v2 :" << n_2 << std::endl;
    std::cout << "valence v3 :" << n_3 << std::endl;

    W = A.inverse() * B;
}   

void WT_Lifting_6(
    const Eigen::Vector3d v_i,
    const double omega_i,
    const Eigen::Vector3d v,
    Eigen::Vector3d& v_i_prime
) {
    v_i_prime = v_i - omega_i * v;
}

void WT_Lifting_6_inverse(
    const Eigen::Vector3d v_i_prime,
    const double omega_i,
    const Eigen::Vector3d v,
    Eigen::Vector3d& v_i
) {
    v_i = v_i_prime + omega_i * v;
}