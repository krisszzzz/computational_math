#include <boost/hana.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/type_index.hpp>
#include <cmath>
#include <matplot/matplot.h>

namespace conv_eq {
struct Ranges {
  int x_steps;
  int t_steps;
  double x_max;
  double t_max;
};

namespace solv_scheme {
struct LaxWendroff {
  void operator()(matplot::vector_2d &U, const Ranges &ranges) const {
    double x_step = ranges.x_max / ranges.x_steps;
    double t_step = ranges.t_max / ranges.t_steps;

    for (int k = 0; k < ranges.t_steps; k++) {
      int m = 1;
      for (; m < ranges.x_steps; m++) {
        // Lax-Wendroff difference scheme
        U[k + 1][m] = U[k][m] -
                      (t_step / (2 * x_step)) * (U[k][m + 1] - U[k][m - 1]) +
                      (0.5 * t_step * t_step / (x_step * x_step)) *
                          (U[k][m + 1] - 2 * U[k][m] + U[k][m - 1]);
      }
      // 3-point difference scheme for boundary
      U[k + 1][m] = U[k][m] - (t_step / x_step) * (U[k][m] - U[k][m - 1]);
    }
  }
};

struct Cross {
  void operator()(matplot::vector_2d &U, const Ranges &ranges) const {
    double x_step = ranges.x_max / ranges.x_steps;
    double t_step = ranges.t_max / ranges.t_steps;
    for (int m = 1; m < ranges.x_steps + 1; m++) {
      // 3-point difference scheme at start (k = 0 => t = 0)
      U[0 + 1][m] = U[0][m] - (t_step / x_step) * (U[0][m] - U[0][m - 1]);
    }

    for (int k = 1; k < ranges.t_steps; k++) {
      int m = 1;
      for (; m < ranges.x_steps; m++) {
        U[k + 1][m] =
            U[k - 1][m] - (t_step / x_step) * (U[k][m + 1] - U[k][m - 1]);
      }
      // 3-point difference scheme for boundary
      U[k + 1][m] = U[k][m] - (t_step / x_step) * (U[k][m] - U[k][m - 1]);
    }
  }
};

struct ThreePoint {
  void operator()(matplot::vector_2d &U, const Ranges &ranges) const {
    double x_step = ranges.x_max / ranges.x_steps;
    double t_step = ranges.t_max / ranges.t_steps;

    for (int k = 0; k < ranges.t_steps; k++) {
      for (int m = 1; m < ranges.x_steps + 1; m++) {
        // 3-point difference scheme
        U[k + 1][m] = U[k][m] - (t_step / x_step) * (U[k][m] - U[k][m - 1]);
      }
    }
  }
};
} // namespace solv_scheme

class SolverBase {
public:
  using InitFunc = std::function<double(double)>;
  // getters
  template <typename Self> auto &&getRanges(this Self &&self) {
    return std::forward<Self>(self).ranges_;
  }

  template <typename Self> auto &&getXCond(this Self &&self) {
    return std::forward<Self>(self).x_cond_;
  }

  template <typename Self> auto &&getTCond(this Self &&self) {
    return std::forward<Self>(self).t_cond_;
  }

  template <typename Self, typename Callable = solv_scheme::ThreePoint>
  matplot::vector_2d solveEquation(this Self &&self,
                                   Callable &&scheme = Callable{}) {
    return std::forward<Self>(self).solveEquationImpl(
        std::forward<Callable>(scheme));
  }

  template <typename Self> auto meshGrid(this Self &&self) {
    return matplot::meshgrid(matplot::linspace(0, self.getRanges().x_max,
                                               self.getRanges().x_steps + 1),
                             matplot::linspace(0, self.getRanges().t_max,
                                               self.getRanges().t_steps + 1));
  }

protected:
  void setInitConds(const InitFunc &x_func, const InitFunc &t_func) {
    double x_step = getRanges().x_max / getRanges().x_steps;
    double t_step = getRanges().t_max / getRanges().t_steps;
    for (int i = 0; i < getRanges().x_steps + 1; i++) {
      getXCond()[i] = x_func(x_step * i);
    }
    for (int i = 0; i < getRanges().t_steps + 1; i++) {
      getTCond()[i] = t_func(t_step * i);
    }
  }

private:
  Ranges ranges_;
  matplot::vector_1d x_cond_;
  matplot::vector_1d t_cond_;
};

// sequantional solver
class SolverSeq : public SolverBase {
public:
  friend class SolverBase;

  SolverSeq(Ranges &ranges, const InitFunc &x_func, const InitFunc &t_func) {
    getRanges() = ranges;
    getXCond().resize(ranges.x_steps + 1);
    getTCond().resize(ranges.t_steps + 1);
    setInitConds(std::forward<const InitFunc>(x_func),
                 std::forward<const InitFunc>(t_func));
  }

private:
  template <typename Callable>
  matplot::vector_2d solveEquationImpl(Callable &&scheme) const {
    matplot::vector_2d U(getRanges().t_steps + 1,
                         matplot::vector_1d(getRanges().x_steps + 1));
    // set init conditions
    for (int i = 0; i < getRanges().x_steps + 1; i++) {
      U[0][i] = getXCond()[i];
    }
    for (int i = 0; i < getRanges().t_steps + 1; i++) {
      U[i][0] = getTCond()[i];
    }
    scheme(U, getRanges());

    return U;
  }
};
} // namespace conv_eq

double calcMaxDiscrepancy(const matplot::vector_2d &analytical_U,
                          const matplot::vector_2d &U) {
  double max_discrepancy = 0.0;

  for (int i = 0; i < analytical_U.size(); i++) {
    for (int j = 0; j < analytical_U[i].size(); j++) {
      double curr_dicrepancy = std::abs(analytical_U[i][j] - U[i][j]);
      max_discrepancy = std::max(curr_dicrepancy, max_discrepancy);
    }
  }

  return max_discrepancy;
}

int main(int argc, char *argv[]) {
  using namespace conv_eq;
  // Solve equation Ut + Ux = 0
  // Boundary U(x, 0) = x
  //          U(0, t) = -t
  // Analytical solution U(x, t) = x - t

  Ranges ranges{
      .x_steps = 40,
      .t_steps = 80,

      .x_max = 3,
      .t_max = 3,
  };
  SolverSeq solver_seq{ranges, [](double x) { return x; },
                       [](double t) { return -t; }};
  double x_step = ranges.x_max / ranges.x_steps;
  double t_step = ranges.t_max / ranges.t_steps;
  matplot::vector_2d analytical_U;

  for (double t = 0; t <= ranges.t_max; t += t_step) {
    matplot::vector_1d U_row;
    for (double x = 0; x <= ranges.x_max; x += x_step) {
      U_row.push_back(x - t);
    }
    analytical_U.push_back(U_row);
  }
  auto tp =
      boost::hana::make_tuple(solv_scheme::LaxWendroff{},
                              solv_scheme::ThreePoint{}, solv_scheme::Cross{});
  boost::hana::for_each(tp, [&solver_seq, &analytical_U](auto &&scheme) {
    constexpr int test_count = 1;
    auto [X, T] = solver_seq.meshGrid();
    matplot::vector_2d U;

    for (int i = 0; i < test_count; i++) {
      U = solver_seq.solveEquation(scheme);
    }

    double discrepancy = calcMaxDiscrepancy(analytical_U, U);
    std::cout << std::format("Implemenation: {:15}, discrepancy = {:.5}",
                             boost::typeindex::type_id<decltype(scheme)>()
                                 .pretty_name(),
                             discrepancy)
              << "\n";

    matplot::surf(X, T, U);
    matplot::show();
  });
}
