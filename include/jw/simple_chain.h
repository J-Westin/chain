#ifndef SIMPLE_CHAIN_H
#define SIMPLE_CHAIN_H



typedef Eigen::Array <double,   N,   1> ArrayNd;
typedef Eigen::Array <double, 2*N,   1> Array2Nd;

typedef Eigen::Matrix<double,   N,   1> VectorNd;
typedef Eigen::Matrix<double, 2*N,   1> Vector2Nd;

typedef Eigen::Matrix<double,   N,   N> MatrixNd;
typedef Eigen::Matrix<double, 2*N, 2*N> Matrix2Nd;

ArrayNd curvy_thetas(){
    ArrayNd thetas(N);

    for ( unsigned int i = 1; i <= N; i++ ){
        thetas[i-1] = i*M_PI/(2.*N);
    }

    return thetas;
}

Matrix2Nd chain_A(ArrayNd thetas){
    ArrayNd s = thetas.sin();
    ArrayNd c = thetas.cos();

    Matrix2Nd A = Matrix2Nd::Zero();

    for (unsigned int i = 0; i < N; i++){
        for (unsigned int j = 0; j < N; j++){
            if (i >= j){
                A(i  , j) =  c(j);
                A(N+i, j) = -s(j);
            }
            if (i == j){
                A(i  , N+j) = s(j);
                A(N+i, N+j) = c(j);
            }
            if (j == i+1){
                A(i  , N+j) = -s(j);
                A(N+i, N+j) = -c(j);
            }
        }
    }

    return A;
}

Vector2Nd chain_b(ArrayNd thetas, ArrayNd omegas){
    ArrayNd s = thetas.sin();
    ArrayNd c = thetas.cos();
    ArrayNd w = omegas.square();

    Vector2Nd b;

    double bx, by;

    for (unsigned int i = 0; i < N; i++){
        bx = 0.;
        by = gamma;

        for (unsigned int j = 0; j <= i; j++){
            bx += w(j)*s(j);
            by += w(j)*c(j);
        }

        b(i)   = bx;
        b(N+i) = by;
    }

    return b;
}

ArrayNd solve_for_alpha(ArrayNd thetas, ArrayNd omegas){
    Matrix2Nd A = chain_A(thetas);
    Vector2Nd b = chain_b(thetas, omegas);

    Vector2Nd x = Eigen::Inverse<Matrix2Nd>(A)*b;
    ArrayNd alphas;

    for (unsigned int i = 0; i < N; i++){
        alphas[i] = x[i];
    }

    return alphas;
}

double modulate_theta(double theta){
    int full_thetas = theta/(2.*M_PI);
    return theta - full_thetas*(2.*M_PI);
}

double clip_omega(double omega, double omega_max=1.){
    if (std::abs(omega) > omega_max){
        if (omega < 0.){
            omega = -omega_max;
        }
        else{
            omega = omega_max;
        }
    }
    return omega;
}

void euler_integrate(ArrayNd& thetas, ArrayNd& omegas, const double dt){
    ArrayNd alphas = solve_for_alpha(thetas, omegas);

    ArrayNd d_theta = dt*omegas;
    ArrayNd d_omega = dt*alphas;

    thetas += d_theta;
    omegas += d_omega;

    for (unsigned int i = 0; i < N; i++){
        thetas[i] = modulate_theta(thetas[i]);
        omegas[i] = clip_omega(omegas[i]);
    }

    return;
}

#endif // SIMPLE_CHAIN_H
