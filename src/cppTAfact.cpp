#include <iostream>
#include <math.h>
#include <limits>
#include <algorithm>
#include <functional>
#include <tuple>
#include <vector>
#include <type_traits>

#include <Eigen/Dense>
#include <Eigen/Cholesky>

/* R-C++ interface with Eigen support */
#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace RcppEigen;

/*
 * Rcpp declarations
 */
// [[Rcpp::depends(RcppEigen)]]
//
// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]


/* Signal handling */
#include <signal.h>
#include <unistd.h>

/* to make Eigen thread-safe */
#include <Eigen/Core>

using Eigen::Map;
using Eigen::Dynamic;
using Eigen::Infinity;
using RMatrixIn  = Map<Eigen::MatrixXd>;
using RMatrixOut = Eigen::MatrixXd;

/* Aliases */
using Double = double;

/* Signal handing */
bool gotSignal = false;
void setGotSignal(int signum) {
    gotSignal = true;
}

/* Binary operator to get projected gradient
 * while optimizing wrt T
 */
class ProjGradT {
    Double eps;

public:
    ProjGradT(Double tol) : eps(tol)
    {}

    inline Double operator()(const Double& g, const Double& t) const {
        if (t <= eps) {
            return std::min(g, 0.0);
        }
        else if (eps < t && t < 1 - eps) {
            return g;
        }
        else {
            return std::max(g, 0.0);
        }
    }
};

template <typename MatrixBig, int DIM = 16, typename Scalar = Double>
class ProbSimplexProjector {
public:
    using Matrix = Eigen::Matrix<Scalar, DIM, Dynamic>;

private:
    const Matrix mTtD;
    const Matrix mTtT;
    double tol;
    int itersMax;

    int niter;
    double optCond;

    int r;
    int n;

public:
    ProbSimplexProjector(const MatrixBig& Dt, const Matrix& Tt, double tol, int itersMax)
        : mTtD(Tt * Dt.transpose()), mTtT(Tt * Tt.transpose()), tol(tol), itersMax(itersMax),
        r(Tt.rows()), n(Dt.rows())
    {}

    void solve(Matrix& mA) {
        /* init */
        niter = 1;
        optCond = 1e+10;

        double cL = mTtT.operatorNorm() + tol;
        double lrA = 1.0 / cL;

        Matrix mAy = mA;
        Matrix mAnext = Matrix::Zero(r, n);
        Matrix gradA  = Matrix::Zero(r, n);
        double tcurr = 1.0, tnext = 1.0;

        while (niter <= itersMax && optCond > tol) {
            evalGrad(mAy, gradA);
            mAnext = mAy - lrA * gradA;
            colwiseProjProbSplx(mAnext);

            /* Check for restart. Gradient-mapping based test */
            if ((cL * (mAy - mAnext)).cwiseProduct(mAnext - mA).sum() > 0) {
                /* Restart */
                mAy = mAnext;
            }
            else {
                mAy = mAnext + (tcurr - 1.0) / tnext * (mAnext - mA);
            }

            tnext = 0.5 * (1 + std::sqrt(1 + 4 * tcurr * tcurr));

            /* Stopping criteria */
            ++niter;
            optCond = (mAnext - mA).norm();

            mA = mAnext;
            tcurr = tnext;
        }
    }

    inline int getNumIters() const {
        return niter;
    }

    inline double getOptCond() const {
        return optCond;
    }

private:
    void evalGrad(const Matrix& A, Matrix& grad) {
        grad = mTtT * A - mTtD;
    }

    void colwiseProjProbSplx(Matrix& mA) {
        size_t n = mA.cols();
        size_t r = mA.rows();
        Matrix mAcopy = mA;

        #pragma omp parallel for
        for (size_t colN = 0; colN < n; ++colN) {
            auto s = mAcopy.col(colN);
            std::sort(s.data(), s.data() + s.size(),
                    std::greater<Double>());
            bool bget = false;
            double tmpsum = 0.0, tmax;
            for (size_t ii = 0; ii < r - 1; ++ii) {
                tmpsum += s(ii);
                tmax = (tmpsum - 1.0) / (ii + 1);
                if (tmax >= s(ii + 1)) {
                    bget = true;
                    break;
                }
            }

            if (!bget) {
                tmax = (tmpsum + s(r - 1) - 1) / r;
            }

            for (size_t jj = 0; jj < r; ++jj) {
                mA(jj, colN) = std::max(mA(jj, colN) - tmax, 0.0);
            }
        }
    }
};

template <int DIM = 16, typename Scalar = Double>
class QPBoxSolverSmallDims {
public:
    using Matrix      = Eigen::Matrix<Scalar, DIM, DIM>;
    using Vector      = Eigen::Matrix<Scalar, DIM, 1>;
    using VectorIdx   = Eigen::Matrix<int, DIM, 1>;
    using StateMatrix = Eigen::Matrix<int, DIM, 3>;

private:
    const Matrix AAt;
    const Vector b;
    double tol;
    int itersMax;
    int r;

    int niter;
    double optCond;

    const Scalar slackEps = 1e-15;

    Scalar L = 1e+05;

public:
    QPBoxSolverSmallDims(const Matrix& AAt, const Vector& b, double tol, int itersMax)
        : AAt(AAt), b(b), tol(tol), itersMax(itersMax), r(AAt.cols())
    {}

    enum Method {newton = 0, coord_descent, fista, exact_any_rank, exact_rank_2};

    void solve(Vector& tinit, int method = newton) {
        tinit = tinit.cwiseMax(0.0).cwiseMin(1.0);
        niter = 1;
        optCond = 1e+10;

        switch(method) {
            case newton:
                solveNewton(tinit);
                break;
            case coord_descent:
                solveCoordDescent(tinit);
                break;
            case fista:
                solveFista(tinit);
                break;
            case exact_any_rank:
                solveExactAnyRank(tinit);
                break;
            case exact_rank_2:
                if (r == 2) {
                    solveExactRank2(tinit);
                }
                else {
                    solveExactAnyRank(tinit);
                }
                break;
            default:
                solveNewton(tinit);
        }
    }

    inline int getNumIters() const {
        return niter;
    }

    inline double getOptCond() const {
        return optCond;
    }

    void setLipschitzConstant(Scalar mL) {
        L = mL;
    }

private:
    enum varState {box = 0, zero, one};

    void solveExactAnyRank(Vector& t) {
        VectorIdx stateVec = VectorIdx::Zero();
        StateMatrix stateOrder = StateMatrix::Zero();
        /* predicting states given initial value */
        for (int i = 0; i < r; ++i) {
            if (t(i) <= slackEps) {
                stateOrder(i, 0) = zero;
                stateOrder(i, 1) = box;
                stateOrder(i, 2) = one;
            }
            else if (t(i) >= 1.0 - slackEps) {
                stateOrder(i, 0) = one;
                stateOrder(i, 1) = box;
                stateOrder(i, 2) = zero;
            }
            else {
                if (t(i) < 0.5) {
                    stateOrder(i, 0) = box;
                    stateOrder(i, 1) = zero;
                    stateOrder(i, 2) = one;
                }
                else {
                    stateOrder(i, 0) = box;
                    stateOrder(i, 1) = one;
                    stateOrder(i, 2) = zero;
                }
            }
        }

        exhaustiveKKT(t, stateVec, stateOrder, 0);
    }

    bool exhaustiveKKT(Vector& t, VectorIdx& stateVec,
            StateMatrix& stateOrder, int currVar) {
        if (currVar == r) {
            /* Maps state to index */
            VectorIdx stateIdx[] = {VectorIdx::Zero(), VectorIdx::Zero(), VectorIdx::Zero()};
            int numState[] = {0, 0, 0};

            Vector colSumStateOne = Vector::Zero();
            for (int i = 0; i < r; ++i) {
                /* which variables are in state i? And how many of them? */
                stateIdx[stateVec(i)](numState[stateVec(i)]++) = i;
                if (one == stateVec(i)) {
                    colSumStateOne += AAt.col(i);
                }
            }

            /* now we prepare a system of linear equations to figure out
             * what would the values for box variables be like */
            Vector tbox = Vector::Zero();
            if (numState[box] > 0) {
                Matrix AAtBox = Matrix::Zero();
                Vector bBox   = Vector::Zero();
                for (int j = 0; j < numState[box]; ++j) {
                    for (int i = 0; i < numState[box]; ++i) {
                        AAtBox(i, j) = AAt(stateIdx[box](i), stateIdx[box](j));
                    }
                    bBox(j) = colSumStateOne(stateIdx[box](j)) - b(stateIdx[box](j));
                }
                /* solve the system */
                tbox = -AAtBox.ldlt().solve(bBox);
                /*
                * check the feasibility for variables
                * corresponding to the box contraints
                */
                for (int i = 0; i < numState[box]; ++i) {
                    //TODO: relax to >= and >=
                    if (tbox(i) <= 0.0 || tbox(i) >= 1.0) {
                        return false;
                    }
                }
            }

            /*
             * find A * t product
             */
            Vector prod = colSumStateOne;
            for (int i = 0; i < numState[box]; ++i) {
                prod += AAt.col(stateIdx[box](i)) * tbox(i);
            }

            /*
             * now that we have A * t product and we know all the values of
             * current solution, it is time to check KKT for the 0-1 components
             */

            /*
             * KKT for 0
             */
            for (int i = 0; i < numState[zero]; ++i) {
                if (prod(stateIdx[zero](i)) <= b(stateIdx[zero](i))) {
                    return false;
                }
            }

            /*
             * KKT for 1
             */
            for (int i = 0; i < numState[one]; ++i) {
                if (prod(stateIdx[one](i)) >= b(stateIdx[one](i))) {
                    return false;
                }
            }

            /*
             * Forming output
             */
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < numState[i]; ++j) {
                    switch (i) {
                        case box:
                            t(stateIdx[i](j)) = tbox(j);
                            break;
                        case zero:
                            t(stateIdx[i](j)) = 0.0;
                            break;
                        case one:
                            t(stateIdx[i](j)) = 1.0;
                            break;
                    }
                }
            }

            return true;
        }

        for (int i = 0; i < 3; ++i) {
            stateVec(currVar) = stateOrder(currVar, i);
            if (exhaustiveKKT(t, stateVec, stateOrder, currVar + 1)) {
                return true;
            }
        }

        return false;
    }

    void solveExactRank2(Vector& t) {
        Double t0;
        Double t1;
        if (b(0) <= 0 && b(1) <= 0) {
            t0 = t1 = 0.0;
        }
        else if (b(0) >= AAt(0, 0) + AAt(0, 1)
                && b(1) >= AAt(1, 0) + AAt(1, 1)) {
            t0 = t1 = 1.0;
        }
        else if (b(0) * AAt(1, 1) > b(1) * AAt(0, 1)
                && b(1) * AAt(0, 0) > b(0) * AAt(0, 1)
                && b(0) * AAt(1, 1) < AAt(0, 0) * AAt(1, 1)
                                    - AAt(0, 1) * AAt(0, 1)
                                    + AAt(0, 1) * b(1)
                && b(1) * AAt(0, 0) < AAt(0, 0) * AAt(1, 1)
                                    - AAt(0, 1) * AAt(0, 1)
                                    + AAt(0, 1) * b(0)) {
            Double invd = 1.0 / (AAt(0, 0) * AAt(1, 1) - AAt(0, 1) * AAt(0, 1));
            t0 = (AAt(1, 1) * b(0) - AAt(0, 1) * b(1)) * invd;
            t1 = (AAt(0, 0) * b(1) - AAt(1, 0) * b(0)) * invd;
        }
        else if (b(1) > 0 && b(1) < AAt(1, 1)
                && b(0) * AAt(1, 1) <= b(1) * AAt(0, 1)) {
            t0 = 0.0;
            t1 = b(1) / AAt(1, 1);
        }
        else if (b(0) > 0 && b(0) < AAt(0, 0)
                && b(1) * AAt(0, 0) <= b(0) * AAt(0, 1)) {
            t0 = b(0) / AAt(0, 0);
            t1 = 0.0;
        }
        else if (b(1) * AAt(0, 0) >= AAt(0, 0) * AAt(1, 1)
                                    - AAt(0, 1) * AAt(0, 1)
                                    + AAt(0, 1) * b(0)
                && b(0) * AAt(1, 1) > AAt(1, 0) * AAt(1, 1)
                && b(0) * AAt(1, 1) < AAt(0, 1) * AAt(1, 1)
                                    + AAt(0, 0) * AAt(1, 1)) {
            t0 = (b(0) - AAt(0, 1)) / AAt(0, 0);
            t1 = 1.0;
        }
        else if (b(0) * AAt(1, 1) >= AAt(0, 0) * AAt(1, 1)
                                    - AAt(0, 1) * AAt(0, 1)
                                    + AAt(0, 1) * b(1)
                && b(1) * AAt(0, 0) > AAt(0, 1) * AAt(0, 0)
                && b(1) * AAt(0, 0) < AAt(0, 1) * AAt(0, 0)
                + AAt(0, 0) * AAt(1, 1)) {
            t0 = 1.0;
            t1 = (b(1) - AAt(0, 1)) / AAt(1, 1);
        }
        else if (b(0) * AAt(1, 1) <= AAt(1, 1) * AAt(0, 1)
                && b(1) * AAt(0, 0) >= AAt(0, 0) * AAt(1, 1)) {
            t0 = 0.0;
            t1 = 1.0;
        }
        else {
            t0 = 1.0;
            t1 = 0.0;
        }

        t(0) = t0;
        t(1) = t1;
    }

    void solveFista(Vector& t) {
        Vector grad;
        Vector ty = t;
        Vector tx;
        Scalar tcurr = 1.0, tnext = 1.0;
        Scalar step = 1.0 / L;

        while (niter <= itersMax && optCond > tol) {
            evalGrad(ty, grad);
            tx = (ty - step * grad).cwiseMax(0.0).cwiseMin(1.0);
            tnext = 0.5 * (1.0 + std::sqrt(1.0 + 4.0 * tcurr * tcurr));
            ty = tx + (tcurr - 1.0) / tnext * (tx - t);

            t = tx;
            tcurr = tnext;

            optCond = grad.binaryExpr(t, ProjGradT(slackEps)).norm();

            /* finish this iteration */
            ++niter;
        }
    }

    void solveCoordDescent(Vector& t) {
        Vector grad;
        Scalar prod;
        optCond = tol + 1.0;
        while (niter <= itersMax && optCond > tol) {
            grad = -b;
            for (int i = 0; i < r; ++i) {
                prod = AAt.col(i).dot(t) - b(i);
                t(i) -= prod / AAt(i, i);
                t(i) = std::max(0.0, t(i));
                t(i) = std::min(1.0, t(i));

                grad += AAt.col(i) * t(i);
            }

            /* Evaluate optimality condition */
            optCond = grad.binaryExpr(t, ProjGradT(slackEps)).norm();

            /* finish this iteration */
            ++niter;
        }
    }

    void solveNewton(Vector& t) {
        VectorIdx fidxv;
        Vector grad, descent;
        evalGrad(t, grad);

        optCond = grad.binaryExpr(t, ProjGradT(slackEps)).norm();
        while (niter <= itersMax && optCond > tol) {
            /* find 'free' and 'restricted' variables */
            int freeVarsNum = 0;
            for (int i = 0; i < r; ++i) {
                /* specifying descent for 'restricted' variables */
                if (t(i) <= slackEps && grad(i) > 0) {
                    descent(i) = 0.0 - t(i);
                }
                else if (t(i) >= 1.0 - slackEps && grad(i) < 0) {
                    descent(i) = 1.0 - t(i);
                }
                else {
                    fidxv(freeVarsNum++) = i;
                }
            }

            /*
            * IMPORTANT
            * if all variables are restricted,
            * then it makes no sense to proceed.
            */
            if (!freeVarsNum) break;

            /* compute the corresponding submatrix of the Hessian */
            Matrix AAtFree = Matrix::Zero();
            for (int j = 0; j < freeVarsNum; ++j) {
                for (int i = 0; i < freeVarsNum; ++i) {
                    AAtFree(i, j) = AAt(fidxv(i), fidxv(j));
                }
            }

            Vector gradFree = Vector::Zero();
            for (int i = 0; i < freeVarsNum; ++i) {
                gradFree(i) = grad(fidxv(i));
            }

            /* Descent for 'free' variables */
            Vector descentFree = -AAtFree.ldlt().solve(gradFree);
            for (int i = 0; i < freeVarsNum; ++i) {
                descent(fidxv(i)) = descentFree(i);
            }

            applyDescent(t, grad, descent);
            evalGrad(t, grad);

            /* Evaluate optimality condition */
            optCond = std::min(grad.binaryExpr(t, ProjGradT(slackEps)).norm(),
                    1e+03 * descent.norm());

            /* finish this iteration */
            ++niter;
        }
    }

    inline Scalar objF(const Vector& t) const {
        return 0.5 * t.dot(AAt * t) - t.dot(b);
    }

    inline void evalGrad(const Vector& t, Vector& g) {
        g = AAt * t - b;
    }

    Scalar applyDescent(Vector& t,
            const Vector& gradCurr,
            const Vector& descent,
            Scalar alpha = 0.1,
            Scalar beta  = 0.5) {
        Scalar objFtcurr = objF(t);
        Vector tcurr = t;
        Vector& tnew = t;

        Scalar step = 1.0;
        Scalar objFtnew;
        while (true) {
            tnew = (tcurr + step * descent).cwiseMin(1.0).cwiseMax(0.0);
            objFtnew = objF(tnew);

            if (objFtnew > objFtcurr + alpha * step * gradCurr.dot(tnew - tcurr)) {
                step *= beta;
            }
            else {
                break;
            }
        }

        return step;
    }
};

struct SolverSuppOutput {
    int niters;
    double objF;
    double rmse;
};

template <int DIM = -1>
void applySolver(const RMatrixIn& Dt, const RMatrixIn& mTtinit, const RMatrixIn& mAinit,
        double lambda, int itersMax, double tol, double tolA, double tolT,
        RMatrixOut& mTtout, RMatrixOut& mAout, SolverSuppOutput& supp) {
    using MatrixDD = Eigen::Matrix<Double, DIM, DIM>;
    using VectorDD = Eigen::Matrix<Double, DIM, 1>;
    using MatrixDX = Eigen::Matrix<Double, DIM, Dynamic>;

    size_t r = mAinit.rows();
    size_t n = Dt.rows();
    size_t m = Dt.cols();

    /* Convert to Eigen's data types */
    MatrixDX Tt = mTtinit.template cast<Double>();
    MatrixDX A  = mAinit.template cast<Double>();

    /* Time-savers */
    auto onesrm = MatrixDX::Ones(r, m);

    //TODO: make it a parameter!!!
    int innerItersMax = 500;

    int niter = 1;
    double optCond = 1e+10;
    int method = QPBoxSolverSmallDims<DIM>::Method::newton;

    MatrixDX Ttprev;
    MatrixDX Aprev;

    while (niter <= itersMax && optCond > tol) {
        Ttprev = Tt;
        Aprev  = A;

        /*
        * Optimization wrt A {
        */
        ProbSimplexProjector<RMatrixIn, DIM> probSmplxProjector(Dt, Tt, tolA, innerItersMax);
        probSmplxProjector.solve(A);

        /*
        * }
        */

        /*
        * Optimization wrt T {
        */
        MatrixDD AAt = A * A.transpose();
        MatrixDX B = A * Dt - lambda * (onesrm - 2 * Ttprev);

        #pragma omp parallel for schedule(runtime)
        for (int i = 0; i < m; ++i) {
            VectorDD t = Tt.col(i);
            VectorDD b = B.col(i);
            QPBoxSolverSmallDims<DIM> solver(AAt, b, tolT, innerItersMax);

            if (2 == r) {
                method = QPBoxSolverSmallDims<DIM>::Method::exact_rank_2;
            }
            else if (14 < r) {
                method = QPBoxSolverSmallDims<DIM>::Method::coord_descent;
            }

            solver.solve(t, method);

            Tt.col(i) = t;
        }
        /*
        * }
        */

        ++niter;
        double dA = (Aprev - A).norm() / std::sqrt(r * n);
        double dT = (Ttprev - Tt).norm() / std::sqrt(r * m);
        optCond = std::sqrt(dA * dA + dT * dT);
    }

    /*
     * Forming output
     */
    mTtout  = Tt;
    mAout   = A;
    supp.niters = niter - 1;
    supp.rmse   = 0.5 * (Dt - A.transpose() * Tt).squaredNorm();
    supp.objF   = supp.rmse + lambda * (Tt.sum() - Tt.squaredNorm());
    supp.rmse  /= m;
    supp.rmse  /= n;
}

/* Some Voodoo magic to eliminate
 * long switches for different dimensions */
template <int ...> struct DimList {};

/* border case */
void solve(int d, const RMatrixIn& mDt, const RMatrixIn& mTtinit, const RMatrixIn& mAinit,
        double lambda, int itersMax, double tol, double tolA, double tolT,
        RMatrixOut& mTtout, RMatrixOut& mAout, SolverSuppOutput& supp,
        DimList<>) {
}

template <int DIM, int ...DIMS>
void solve(int d, const RMatrixIn& mDt, const RMatrixIn& mTtinit, const RMatrixIn& mAinit,
        double lambda, int itersMax, double tol, double tolA, double tolT,
        RMatrixOut& mTtout, RMatrixOut& mAout, SolverSuppOutput& supp,
        DimList<DIM, DIMS...>) {
    if (DIM != d) {
        return solve(d, mDt, mTtinit, mAinit, lambda,
                itersMax, tol, tolA, tolT,
                mTtout, mAout, supp,
                DimList<DIMS...>());
    }

    applySolver<DIM>(mDt, mTtinit, mAinit, lambda,
            itersMax, tol, tolA, tolT,
            mTtout, mAout, supp);
}

template <int ...DIMS>
void solve(int d, const RMatrixIn& mDt, const RMatrixIn& mTtinit, const RMatrixIn& mAinit,
        double lambda, int itersMax, double tol, double tolA, double tolT,
        RMatrixOut& mTtout, RMatrixOut& mAout, SolverSuppOutput& supp) {
        solve(d, mDt, mTtinit, mAinit, lambda,
                itersMax, tol, tolA, tolT,
                mTtout, mAout, supp,
                DimList<DIMS...>());
}

// [[Rcpp::export]]
RcppExport SEXP cppTAfact(SEXP mDtSEXP, SEXP mTtinitSEXP, SEXP mAinitSEXP,
        double lambda = 0.0, int itersMax = 1000,
        double tol = 1e-8, double tolA = 1e-7, double tolT = 1e-7) {
    /* Prepare Eigen for multithreading */
    Eigen::initParallel();
    Eigen::setNbThreads(1);

    /*
     * We have to set global variables after each call.
     */
    gotSignal = false;
    signal(SIGINT, setGotSignal);
    signal(SIGTERM, setGotSignal);
    signal(SIGKILL, setGotSignal);

    RMatrixIn mDt(as<RMatrixIn>(mDtSEXP));
    RMatrixIn mTtinit(as<RMatrixIn>(mTtinitSEXP));
    RMatrixIn mAinit(as<RMatrixIn>(mAinitSEXP));

    /* Dimensionality of a problem */
    const size_t r = mAinit.rows();

    RMatrixOut mTtout, mAout;
    SolverSuppOutput supp;
    solve<2, 3, 4, 5,
          6, 7, 8, 9,
          10, 11, 12,
          14, 15, 16>(r, mDt, mTtinit, mAinit, lambda, itersMax,
                  tol, tolA, tolT,
                  mTtout, mAout, supp);

    return wrap(List::create(Named("Tt")    = mTtout,
                             Named("A")     = mAout,
                             Named("niter") = supp.niters,
                             Named("objF")  = supp.objF,
                             Named("rmse")  = supp.rmse));
}
