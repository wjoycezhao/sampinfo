// Generated by rstantools.  Do not edit by hand.

/*
    sampinfo is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    sampinfo is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with sampinfo.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#include <rstan/rstaninc.hpp>
// Code generated by Stan version 2.19.1
#include <stan/model/model_header.hpp>
namespace model_simulate_namespace {
using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::prob_grad;
using namespace stan::math;
static int current_statement_begin__;
stan::io::program_reader prog_reader__() {
    stan::io::program_reader reader;
    reader.add_event(0, 0, "start", "model_simulate");
    reader.add_event(153, 151, "end", "model_simulate");
    return reader;
}
template <typename T1__, typename T2__, typename T3__>
std::vector<Eigen::Matrix<typename boost::math::tools::promote_args<T1__, T2__, T3__>::type, Eigen::Dynamic, 1> >
decision_l(const int& N,
               const T1__& ar_value,
               const Eigen::Matrix<T2__, Eigen::Dynamic, 1>& p_n,
               const Eigen::Matrix<T3__, Eigen::Dynamic, 1>& p_y, std::ostream* pstream__) {
    typedef typename boost::math::tools::promote_args<T1__, T2__, T3__>::type local_scalar_t__;
    typedef local_scalar_t__ fun_return_scalar_t__;
    const static bool propto__ = true;
    (void) propto__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
    int current_statement_begin__ = -1;
    try {
        {
        current_statement_begin__ = 6;
        validate_non_negative_index("prob_d", "N", N);
        validate_non_negative_index("prob_d", "3", 3);
        std::vector<Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1>  > prob_d(3, Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1>(N));
        stan::math::initialize(prob_d, DUMMY_VAR__);
        stan::math::fill(prob_d, DUMMY_VAR__);
        current_statement_begin__ = 8;
        stan::model::assign(prob_d, 
                    stan::model::cons_list(stan::model::index_uni(1), stan::model::nil_index_list()), 
                    add(multiply(ar_value, p_n), multiply((1 - ar_value), add(elt_multiply(p_n, subtract(1, p_y)), elt_multiply(multiply(.5, p_y), p_n)))), 
                    "assigning variable prob_d");
        current_statement_begin__ = 11;
        stan::model::assign(prob_d, 
                    stan::model::cons_list(stan::model::index_uni(2), stan::model::nil_index_list()), 
                    add(multiply(ar_value, subtract(subtract(1, p_y), p_n)), elt_multiply(multiply((1 - ar_value), subtract(1, p_n)), subtract(1, p_y))), 
                    "assigning variable prob_d");
        current_statement_begin__ = 14;
        stan::model::assign(prob_d, 
                    stan::model::cons_list(stan::model::index_uni(3), stan::model::nil_index_list()), 
                    add(multiply(ar_value, p_y), multiply((1 - ar_value), add(elt_multiply(p_y, subtract(1, p_n)), elt_multiply(multiply(.5, p_y), p_n)))), 
                    "assigning variable prob_d");
        current_statement_begin__ = 16;
        return stan::math::promote_scalar<fun_return_scalar_t__>(prob_d);
        }
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}
struct decision_l_functor__ {
    template <typename T1__, typename T2__, typename T3__>
        std::vector<Eigen::Matrix<typename boost::math::tools::promote_args<T1__, T2__, T3__>::type, Eigen::Dynamic, 1> >
    operator()(const int& N,
               const T1__& ar_value,
               const Eigen::Matrix<T2__, Eigen::Dynamic, 1>& p_n,
               const Eigen::Matrix<T3__, Eigen::Dynamic, 1>& p_y, std::ostream* pstream__) const {
        return decision_l(N, ar_value, p_n, p_y, pstream__);
    }
};
int
sign(const int& x, std::ostream* pstream__) {
    typedef double local_scalar_t__;
    typedef int fun_return_scalar_t__;
    const static bool propto__ = true;
    (void) propto__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
    int current_statement_begin__ = -1;
    try {
        {
        current_statement_begin__ = 19;
        int x_sign(0);
        (void) x_sign;  // dummy to suppress unused var warning
        stan::math::fill(x_sign, std::numeric_limits<int>::min());
        current_statement_begin__ = 20;
        if (as_bool(logical_lt(x, 0))) {
            current_statement_begin__ = 20;
            stan::math::assign(x_sign, -(1));
        }
        current_statement_begin__ = 21;
        if (as_bool(logical_eq(x, 0))) {
            current_statement_begin__ = 21;
            stan::math::assign(x_sign, 0);
        }
        current_statement_begin__ = 22;
        if (as_bool(logical_gt(x, 0))) {
            current_statement_begin__ = 22;
            stan::math::assign(x_sign, 1);
        }
        current_statement_begin__ = 23;
        return stan::math::promote_scalar<fun_return_scalar_t__>(x_sign);
        }
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}
struct sign_functor__ {
            int
    operator()(const int& x, std::ostream* pstream__) const {
        return sign(x, pstream__);
    }
};
#include <stan_meta_header.hpp>
class model_simulate : public prob_grad {
private:
        int ar_value;
        int deltaM_value;
        int deltaD_value;
        int C;
        int MC;
        int K;
        std::vector<matrix_d> X_sim;
        int max_tNo_prd;
        std::vector<double> deltaM;
        vector_d alpha;
        vector_d beta;
        std::vector<double> deltaD;
        double threshold;
        double sigma;
public:
    model_simulate(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : prob_grad(0) {
        ctor_body(context__, 0, pstream__);
    }
    model_simulate(stan::io::var_context& context__,
        unsigned int random_seed__,
        std::ostream* pstream__ = 0)
        : prob_grad(0) {
        ctor_body(context__, random_seed__, pstream__);
    }
    void ctor_body(stan::io::var_context& context__,
                   unsigned int random_seed__,
                   std::ostream* pstream__) {
        typedef double local_scalar_t__;
        boost::ecuyer1988 base_rng__ =
          stan::services::util::create_rng(random_seed__, 0);
        (void) base_rng__;  // suppress unused var warning
        current_statement_begin__ = -1;
        static const char* function__ = "model_simulate_namespace::model_simulate";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        try {
            // initialize data block variables from context__
            current_statement_begin__ = 29;
            context__.validate_dims("data initialization", "ar_value", "int", context__.to_vec());
            ar_value = int(0);
            vals_i__ = context__.vals_i("ar_value");
            pos__ = 0;
            ar_value = vals_i__[pos__++];
            current_statement_begin__ = 30;
            context__.validate_dims("data initialization", "deltaM_value", "int", context__.to_vec());
            deltaM_value = int(0);
            vals_i__ = context__.vals_i("deltaM_value");
            pos__ = 0;
            deltaM_value = vals_i__[pos__++];
            current_statement_begin__ = 31;
            context__.validate_dims("data initialization", "deltaD_value", "int", context__.to_vec());
            deltaD_value = int(0);
            vals_i__ = context__.vals_i("deltaD_value");
            pos__ = 0;
            deltaD_value = vals_i__[pos__++];
            current_statement_begin__ = 32;
            context__.validate_dims("data initialization", "C", "int", context__.to_vec());
            C = int(0);
            vals_i__ = context__.vals_i("C");
            pos__ = 0;
            C = vals_i__[pos__++];
            check_greater_or_equal(function__, "C", C, 2);
            current_statement_begin__ = 33;
            context__.validate_dims("data initialization", "MC", "int", context__.to_vec());
            MC = int(0);
            vals_i__ = context__.vals_i("MC");
            pos__ = 0;
            MC = vals_i__[pos__++];
            check_greater_or_equal(function__, "MC", MC, 2);
            current_statement_begin__ = 34;
            context__.validate_dims("data initialization", "K", "int", context__.to_vec());
            K = int(0);
            vals_i__ = context__.vals_i("K");
            pos__ = 0;
            K = vals_i__[pos__++];
            check_greater_or_equal(function__, "K", K, 0);
            current_statement_begin__ = 35;
            validate_non_negative_index("X_sim", "C", C);
            validate_non_negative_index("X_sim", "K", K);
            validate_non_negative_index("X_sim", "C", C);
            context__.validate_dims("data initialization", "X_sim", "matrix_d", context__.to_vec(C,C,K));
            X_sim = std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> >(C, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(C, K));
            vals_r__ = context__.vals_r("X_sim");
            pos__ = 0;
            size_t X_sim_j_2_max__ = K;
            size_t X_sim_j_1_max__ = C;
            size_t X_sim_k_0_max__ = C;
            for (size_t j_2__ = 0; j_2__ < X_sim_j_2_max__; ++j_2__) {
                for (size_t j_1__ = 0; j_1__ < X_sim_j_1_max__; ++j_1__) {
                    for (size_t k_0__ = 0; k_0__ < X_sim_k_0_max__; ++k_0__) {
                        X_sim[k_0__](j_1__, j_2__) = vals_r__[pos__++];
                    }
                }
            }
            current_statement_begin__ = 36;
            context__.validate_dims("data initialization", "max_tNo_prd", "int", context__.to_vec());
            max_tNo_prd = int(0);
            vals_i__ = context__.vals_i("max_tNo_prd");
            pos__ = 0;
            max_tNo_prd = vals_i__[pos__++];
            current_statement_begin__ = 37;
            validate_non_negative_index("deltaM", "2", 2);
            context__.validate_dims("data initialization", "deltaM", "double", context__.to_vec(2));
            deltaM = std::vector<double>(2, double(0));
            vals_r__ = context__.vals_r("deltaM");
            pos__ = 0;
            size_t deltaM_k_0_max__ = 2;
            for (size_t k_0__ = 0; k_0__ < deltaM_k_0_max__; ++k_0__) {
                deltaM[k_0__] = vals_r__[pos__++];
            }
            size_t deltaM_i_0_max__ = 2;
            for (size_t i_0__ = 0; i_0__ < deltaM_i_0_max__; ++i_0__) {
                check_greater_or_equal(function__, "deltaM[i_0__]", deltaM[i_0__], 0);
                check_less_or_equal(function__, "deltaM[i_0__]", deltaM[i_0__], 1);
            }
            current_statement_begin__ = 38;
            validate_non_negative_index("alpha", "C", C);
            context__.validate_dims("data initialization", "alpha", "vector_d", context__.to_vec(C));
            alpha = Eigen::Matrix<double, Eigen::Dynamic, 1>(C);
            vals_r__ = context__.vals_r("alpha");
            pos__ = 0;
            size_t alpha_j_1_max__ = C;
            for (size_t j_1__ = 0; j_1__ < alpha_j_1_max__; ++j_1__) {
                alpha(j_1__) = vals_r__[pos__++];
            }
            current_statement_begin__ = 39;
            validate_non_negative_index("beta", "K", K);
            context__.validate_dims("data initialization", "beta", "vector_d", context__.to_vec(K));
            beta = Eigen::Matrix<double, Eigen::Dynamic, 1>(K);
            vals_r__ = context__.vals_r("beta");
            pos__ = 0;
            size_t beta_j_1_max__ = K;
            for (size_t j_1__ = 0; j_1__ < beta_j_1_max__; ++j_1__) {
                beta(j_1__) = vals_r__[pos__++];
            }
            current_statement_begin__ = 40;
            validate_non_negative_index("deltaD", "2", 2);
            context__.validate_dims("data initialization", "deltaD", "double", context__.to_vec(2));
            deltaD = std::vector<double>(2, double(0));
            vals_r__ = context__.vals_r("deltaD");
            pos__ = 0;
            size_t deltaD_k_0_max__ = 2;
            for (size_t k_0__ = 0; k_0__ < deltaD_k_0_max__; ++k_0__) {
                deltaD[k_0__] = vals_r__[pos__++];
            }
            size_t deltaD_i_0_max__ = 2;
            for (size_t i_0__ = 0; i_0__ < deltaD_i_0_max__; ++i_0__) {
                check_greater_or_equal(function__, "deltaD[i_0__]", deltaD[i_0__], 0);
                check_less_or_equal(function__, "deltaD[i_0__]", deltaD[i_0__], 1);
            }
            current_statement_begin__ = 41;
            context__.validate_dims("data initialization", "threshold", "double", context__.to_vec());
            threshold = double(0);
            vals_r__ = context__.vals_r("threshold");
            pos__ = 0;
            threshold = vals_r__[pos__++];
            check_greater_or_equal(function__, "threshold", threshold, 0);
            current_statement_begin__ = 42;
            context__.validate_dims("data initialization", "sigma", "double", context__.to_vec());
            sigma = double(0);
            vals_r__ = context__.vals_r("sigma");
            pos__ = 0;
            sigma = vals_r__[pos__++];
            check_greater_or_equal(function__, "sigma", sigma, 0);
            // initialize transformed data variables
            // execute transformed data statements
            // validate transformed data
            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    ~model_simulate() { }
    void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__,
                         std::ostream* pstream__) const {
        typedef double local_scalar_t__;
        stan::io::writer<double> writer__(params_r__, params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;
        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
    }
    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream__) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec, pstream__);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }
    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(std::vector<T__>& params_r__,
                 std::vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const {
        typedef T__ local_scalar_t__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // dummy to suppress unused var warning
        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;
        try {
            stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
            // model parameters
            // model body
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
        lp_accum__.add(lp__);
        return lp_accum__.sum();
    } // log_prob()
    template <bool propto, bool jacobian, typename T_>
    T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
               std::ostream* pstream = 0) const {
      std::vector<T_> vec_params_r;
      vec_params_r.reserve(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r(i));
      std::vector<int> vec_params_i;
      return log_prob<propto,jacobian,T_>(vec_params_r, vec_params_i, pstream);
    }
    void get_param_names(std::vector<std::string>& names__) const {
        names__.resize(0);
        names__.push_back("rating_prd");
        names__.push_back("cID_prd");
        names__.push_back("terminate_prd");
        names__.push_back("decision_p_prd1");
        names__.push_back("tNo_prd");
    }
    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dims__.push_back(max_tNo_prd);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(max_tNo_prd);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(max_tNo_prd);
        dims__.push_back(3);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
    }
    template <typename RNG>
    void write_array(RNG& base_rng__,
                     std::vector<double>& params_r__,
                     std::vector<int>& params_i__,
                     std::vector<double>& vars__,
                     bool include_tparams__ = true,
                     bool include_gqs__ = true,
                     std::ostream* pstream__ = 0) const {
        typedef double local_scalar_t__;
        vars__.resize(0);
        stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
        static const char* function__ = "model_simulate_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        if (!include_tparams__ && !include_gqs__) return;
        try {
            if (!include_gqs__ && !include_tparams__) return;
            if (!include_gqs__) return;
            // declare and define generated quantities
            current_statement_begin__ = 52;
            validate_non_negative_index("rating_prd", "max_tNo_prd", max_tNo_prd);
            std::vector<int> rating_prd(max_tNo_prd, int(0));
            stan::math::fill(rating_prd, std::numeric_limits<int>::min());
            current_statement_begin__ = 53;
            validate_non_negative_index("cID_prd", "max_tNo_prd", max_tNo_prd);
            std::vector<int> cID_prd(max_tNo_prd, int(0));
            stan::math::fill(cID_prd, std::numeric_limits<int>::min());
            current_statement_begin__ = 54;
            int terminate_prd;
            (void) terminate_prd;  // dummy to suppress unused var warning
            stan::math::fill(terminate_prd, std::numeric_limits<int>::min());
            current_statement_begin__ = 55;
            validate_non_negative_index("decision_p_prd1", "3", 3);
            validate_non_negative_index("decision_p_prd1", "max_tNo_prd", max_tNo_prd);
            std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1> > decision_p_prd1(max_tNo_prd, Eigen::Matrix<double, Eigen::Dynamic, 1>(3));
            stan::math::initialize(decision_p_prd1, DUMMY_VAR__);
            stan::math::fill(decision_p_prd1, DUMMY_VAR__);
            current_statement_begin__ = 56;
            int tNo_prd;
            (void) tNo_prd;  // dummy to suppress unused var warning
            stan::math::fill(tNo_prd, std::numeric_limits<int>::min());
            // generated quantities statements
            current_statement_begin__ = 59;
            stan::math::assign(cID_prd, rep_array(99, max_tNo_prd));
            current_statement_begin__ = 60;
            stan::math::assign(rating_prd, rep_array(99, max_tNo_prd));
            current_statement_begin__ = 61;
            for (int i = 1; i <= max_tNo_prd; ++i) {
                current_statement_begin__ = 62;
                stan::model::assign(decision_p_prd1, 
                            stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                            to_vector(rep_array(99, 3)), 
                            "assigning variable decision_p_prd1");
            }
            {
            current_statement_begin__ = 66;
            validate_non_negative_index("X_acc_prd", "C", C);
            validate_non_negative_index("X_acc_prd", "K", K);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, Eigen::Dynamic> X_acc_prd(C, K);
            stan::math::initialize(X_acc_prd, DUMMY_VAR__);
            stan::math::fill(X_acc_prd, DUMMY_VAR__);
            current_statement_begin__ = 67;
            validate_non_negative_index("theta_prd", "C", C);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> theta_prd(C);
            stan::math::initialize(theta_prd, DUMMY_VAR__);
            stan::math::fill(theta_prd, DUMMY_VAR__);
            current_statement_begin__ = 68;
            int rating_n_prd(0);
            (void) rating_n_prd;  // dummy to suppress unused var warning
            stan::math::fill(rating_n_prd, std::numeric_limits<int>::min());
            current_statement_begin__ = 69;
            int rating_y_prd(0);
            (void) rating_y_prd;  // dummy to suppress unused var warning
            stan::math::fill(rating_y_prd, std::numeric_limits<int>::min());
            current_statement_begin__ = 70;
            local_scalar_t__ utility_n_prd(DUMMY_VAR__);
            (void) utility_n_prd;  // dummy to suppress unused var warning
            stan::math::initialize(utility_n_prd, DUMMY_VAR__);
            stan::math::fill(utility_n_prd, DUMMY_VAR__);
            current_statement_begin__ = 71;
            local_scalar_t__ utility_y_prd(DUMMY_VAR__);
            (void) utility_y_prd;  // dummy to suppress unused var warning
            stan::math::initialize(utility_y_prd, DUMMY_VAR__);
            stan::math::fill(utility_y_prd, DUMMY_VAR__);
            current_statement_begin__ = 72;
            local_scalar_t__ p_n_prd(DUMMY_VAR__);
            (void) p_n_prd;  // dummy to suppress unused var warning
            stan::math::initialize(p_n_prd, DUMMY_VAR__);
            stan::math::fill(p_n_prd, DUMMY_VAR__);
            current_statement_begin__ = 73;
            local_scalar_t__ p_y_prd(DUMMY_VAR__);
            (void) p_y_prd;  // dummy to suppress unused var warning
            stan::math::initialize(p_y_prd, DUMMY_VAR__);
            stan::math::fill(p_y_prd, DUMMY_VAR__);
            current_statement_begin__ = 74;
            validate_non_negative_index("decision_p_prd", "1", 1);
            validate_non_negative_index("decision_p_prd", "3", 3);
            std::vector<Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1>  > decision_p_prd(3, Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1>(1));
            stan::math::initialize(decision_p_prd, DUMMY_VAR__);
            stan::math::fill(decision_p_prd, DUMMY_VAR__);
            current_statement_begin__ = 76;
            stan::math::assign(tNo_prd, 1);
            current_statement_begin__ = 77;
            stan::math::assign(terminate_prd, 0);
            current_statement_begin__ = 78;
            while (as_bool((primitive_value(logical_eq(terminate_prd, 0)) && primitive_value(logical_lt(tNo_prd, (max_tNo_prd + 1)))))) {
                current_statement_begin__ = 79;
                if (as_bool((primitive_value(logical_gt(K, 0)) && primitive_value(logical_eq(deltaM_value, 8))))) {
                    current_statement_begin__ = 80;
                    if (as_bool(logical_eq(tNo_prd, 1))) {
                        current_statement_begin__ = 81;
                        stan::math::assign(X_acc_prd, to_matrix(rep_array(0, C, K)));
                    } else {
                        current_statement_begin__ = 83;
                        stan::math::assign(X_acc_prd, multiply(X_acc_prd, get_base1(deltaM, 1, "deltaM", 1)));
                        current_statement_begin__ = 84;
                        stan::math::assign(X_acc_prd, add(X_acc_prd, get_base1(X_sim, get_base1(cID_prd, (tNo_prd - 1), "cID_prd", 1), "X_sim", 1)));
                    }
                    current_statement_begin__ = 86;
                    stan::math::assign(theta_prd, multiply(X_acc_prd, beta));
                    current_statement_begin__ = 87;
                    stan::math::assign(theta_prd, add(theta_prd, alpha));
                } else if (as_bool((primitive_value(logical_gt(K, 0)) && primitive_value(logical_eq(deltaM_value, 0))))) {
                    current_statement_begin__ = 89;
                    if (as_bool(logical_eq(tNo_prd, 1))) {
                        current_statement_begin__ = 90;
                        stan::math::assign(X_acc_prd, to_matrix(rep_array(0, C, K)));
                    } else {
                        current_statement_begin__ = 92;
                        stan::math::assign(X_acc_prd, get_base1(X_sim, get_base1(cID_prd, (tNo_prd - 1), "cID_prd", 1), "X_sim", 1));
                    }
                    current_statement_begin__ = 94;
                    stan::math::assign(theta_prd, multiply(X_acc_prd, beta));
                    current_statement_begin__ = 95;
                    stan::math::assign(theta_prd, add(theta_prd, alpha));
                } else if (as_bool((primitive_value(logical_gt(K, 0)) && primitive_value(logical_eq(deltaM_value, 1))))) {
                    current_statement_begin__ = 97;
                    if (as_bool(logical_eq(tNo_prd, 1))) {
                        current_statement_begin__ = 98;
                        stan::math::assign(X_acc_prd, to_matrix(rep_array(0, C, K)));
                    } else {
                        current_statement_begin__ = 100;
                        stan::math::assign(X_acc_prd, add(X_acc_prd, get_base1(X_sim, get_base1(cID_prd, (tNo_prd - 1), "cID_prd", 1), "X_sim", 1)));
                    }
                    current_statement_begin__ = 102;
                    stan::math::assign(theta_prd, multiply(X_acc_prd, beta));
                    current_statement_begin__ = 103;
                    stan::math::assign(theta_prd, add(theta_prd, alpha));
                } else {
                    current_statement_begin__ = 105;
                    stan::math::assign(theta_prd, alpha);
                }
                current_statement_begin__ = 108;
                stan::model::assign(cID_prd, 
                            stan::model::cons_list(stan::model::index_uni(tNo_prd), stan::model::nil_index_list()), 
                            categorical_logit_rng(theta_prd, base_rng__), 
                            "assigning variable cID_prd");
                current_statement_begin__ = 110;
                stan::model::assign(rating_prd, 
                            stan::model::cons_list(stan::model::index_uni(tNo_prd), stan::model::nil_index_list()), 
                            sign((get_base1(cID_prd, tNo_prd, "cID_prd", 1) - MC), pstream__), 
                            "assigning variable rating_prd");
                current_statement_begin__ = 112;
                if (as_bool(logical_eq(ar_value, 0))) {
                    current_statement_begin__ = 113;
                    stan::math::assign(rating_n_prd, std::min(get_base1(rating_prd, tNo_prd, "rating_prd", 1), 0));
                    current_statement_begin__ = 114;
                    stan::math::assign(rating_y_prd, std::max(get_base1(rating_prd, tNo_prd, "rating_prd", 1), 0));
                } else {
                    current_statement_begin__ = 116;
                    stan::math::assign(rating_n_prd, get_base1(rating_prd, tNo_prd, "rating_prd", 1));
                    current_statement_begin__ = 117;
                    stan::math::assign(rating_y_prd, get_base1(rating_prd, tNo_prd, "rating_prd", 1));
                }
                current_statement_begin__ = 120;
                if (as_bool(logical_eq(deltaD_value, 8))) {
                    current_statement_begin__ = 121;
                    if (as_bool(logical_eq(tNo_prd, 1))) {
                        current_statement_begin__ = 122;
                        stan::math::assign(utility_n_prd, rating_n_prd);
                        current_statement_begin__ = 123;
                        stan::math::assign(utility_y_prd, rating_y_prd);
                    } else {
                        current_statement_begin__ = 125;
                        stan::math::assign(utility_n_prd, (utility_n_prd * get_base1(deltaD, 1, "deltaD", 1)));
                        current_statement_begin__ = 126;
                        stan::math::assign(utility_n_prd, (utility_n_prd + rating_n_prd));
                        current_statement_begin__ = 127;
                        stan::math::assign(utility_y_prd, (utility_y_prd * get_base1(deltaD, 1, "deltaD", 1)));
                        current_statement_begin__ = 128;
                        stan::math::assign(utility_y_prd, (utility_y_prd + rating_y_prd));
                    }
                } else if (as_bool(logical_eq(deltaD_value, 0))) {
                    current_statement_begin__ = 131;
                    stan::math::assign(utility_n_prd, rating_n_prd);
                    current_statement_begin__ = 132;
                    stan::math::assign(utility_y_prd, rating_y_prd);
                } else if (as_bool(logical_eq(deltaD_value, 1))) {
                    current_statement_begin__ = 134;
                    if (as_bool(logical_eq(tNo_prd, 1))) {
                        current_statement_begin__ = 135;
                        stan::math::assign(utility_n_prd, rating_n_prd);
                        current_statement_begin__ = 136;
                        stan::math::assign(utility_y_prd, rating_y_prd);
                    } else {
                        current_statement_begin__ = 138;
                        stan::math::assign(utility_n_prd, (utility_n_prd + rating_n_prd));
                        current_statement_begin__ = 139;
                        stan::math::assign(utility_y_prd, (utility_y_prd + rating_y_prd));
                    }
                }
                current_statement_begin__ = 142;
                stan::math::assign(p_n_prd, (1 - normal_cdf(utility_n_prd, -(threshold), sigma)));
                current_statement_begin__ = 143;
                stan::math::assign(p_y_prd, normal_cdf(utility_y_prd, threshold, sigma));
                current_statement_begin__ = 144;
                stan::math::assign(decision_p_prd, decision_l(1, ar_value, to_vector(rep_array(p_n_prd, 1)), to_vector(rep_array(p_y_prd, 1)), pstream__));
                current_statement_begin__ = 145;
                for (int i = 1; i <= 3; ++i) {
                    current_statement_begin__ = 145;
                    stan::model::assign(decision_p_prd1, 
                                stan::model::cons_list(stan::model::index_uni(tNo_prd), stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list())), 
                                get_base1(get_base1(decision_p_prd, i, "decision_p_prd", 1), 1, "decision_p_prd", 2), 
                                "assigning variable decision_p_prd1");
                }
                current_statement_begin__ = 146;
                stan::math::assign(terminate_prd, (categorical_rng(get_base1(decision_p_prd1, tNo_prd, "decision_p_prd1", 1), base_rng__) - 2));
                current_statement_begin__ = 147;
                stan::math::assign(tNo_prd, (tNo_prd + 1));
            }
            }
            current_statement_begin__ = 150;
            stan::math::assign(tNo_prd, (tNo_prd + -(1)));
            // validate, write generated quantities
            current_statement_begin__ = 52;
            size_t rating_prd_k_0_max__ = max_tNo_prd;
            for (size_t k_0__ = 0; k_0__ < rating_prd_k_0_max__; ++k_0__) {
                vars__.push_back(rating_prd[k_0__]);
            }
            current_statement_begin__ = 53;
            size_t cID_prd_k_0_max__ = max_tNo_prd;
            for (size_t k_0__ = 0; k_0__ < cID_prd_k_0_max__; ++k_0__) {
                vars__.push_back(cID_prd[k_0__]);
            }
            current_statement_begin__ = 54;
            vars__.push_back(terminate_prd);
            current_statement_begin__ = 55;
            size_t decision_p_prd1_j_1_max__ = 3;
            size_t decision_p_prd1_k_0_max__ = max_tNo_prd;
            for (size_t j_1__ = 0; j_1__ < decision_p_prd1_j_1_max__; ++j_1__) {
                for (size_t k_0__ = 0; k_0__ < decision_p_prd1_k_0_max__; ++k_0__) {
                    vars__.push_back(decision_p_prd1[k_0__](j_1__));
                }
            }
            current_statement_begin__ = 56;
            vars__.push_back(tNo_prd);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    template <typename RNG>
    void write_array(RNG& base_rng,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                     bool include_tparams = true,
                     bool include_gqs = true,
                     std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<double> vars_vec;
      std::vector<int> params_i_vec;
      write_array(base_rng, params_r_vec, params_i_vec, vars_vec, include_tparams, include_gqs, pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }
    static std::string model_name() {
        return "model_simulate";
    }
    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
        }
        if (!include_gqs__) return;
        size_t rating_prd_k_0_max__ = max_tNo_prd;
        for (size_t k_0__ = 0; k_0__ < rating_prd_k_0_max__; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "rating_prd" << '.' << k_0__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t cID_prd_k_0_max__ = max_tNo_prd;
        for (size_t k_0__ = 0; k_0__ < cID_prd_k_0_max__; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "cID_prd" << '.' << k_0__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        param_name_stream__.str(std::string());
        param_name_stream__ << "terminate_prd";
        param_names__.push_back(param_name_stream__.str());
        size_t decision_p_prd1_j_1_max__ = 3;
        size_t decision_p_prd1_k_0_max__ = max_tNo_prd;
        for (size_t j_1__ = 0; j_1__ < decision_p_prd1_j_1_max__; ++j_1__) {
            for (size_t k_0__ = 0; k_0__ < decision_p_prd1_k_0_max__; ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "decision_p_prd1" << '.' << k_0__ + 1 << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        param_name_stream__.str(std::string());
        param_name_stream__ << "tNo_prd";
        param_names__.push_back(param_name_stream__.str());
    }
    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
        }
        if (!include_gqs__) return;
        size_t rating_prd_k_0_max__ = max_tNo_prd;
        for (size_t k_0__ = 0; k_0__ < rating_prd_k_0_max__; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "rating_prd" << '.' << k_0__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t cID_prd_k_0_max__ = max_tNo_prd;
        for (size_t k_0__ = 0; k_0__ < cID_prd_k_0_max__; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "cID_prd" << '.' << k_0__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        param_name_stream__.str(std::string());
        param_name_stream__ << "terminate_prd";
        param_names__.push_back(param_name_stream__.str());
        size_t decision_p_prd1_j_1_max__ = 3;
        size_t decision_p_prd1_k_0_max__ = max_tNo_prd;
        for (size_t j_1__ = 0; j_1__ < decision_p_prd1_j_1_max__; ++j_1__) {
            for (size_t k_0__ = 0; k_0__ < decision_p_prd1_k_0_max__; ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "decision_p_prd1" << '.' << k_0__ + 1 << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        param_name_stream__.str(std::string());
        param_name_stream__ << "tNo_prd";
        param_names__.push_back(param_name_stream__.str());
    }
}; // model
}  // namespace
typedef model_simulate_namespace::model_simulate stan_model;
#endif
