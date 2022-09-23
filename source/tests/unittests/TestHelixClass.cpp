#include "HelixClass.h"
#include "HelixClass_double.h"

#include "CLHEP/Units/PhysicalConstants.h"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <tuple>

// The speed of light as used inside the Helix class
constexpr auto c_helix = CLHEP::c_light / (CLHEP::m / CLHEP::ps);

// Using this type list and a TEMPLATE_LIST_TEST_CASE we can test both versions
// of the HelixClass with only one set of code. Inside the test case we can
// refer to the current HelixClass type with TestType
using HelixTypes = std::tuple<HelixClass, HelixClass_double>;

TEMPLATE_LIST_TEST_CASE("Initialize_VP", "[helix-init]", HelixTypes) {
  TestType helix;
  using FloatT = typename TestType::float_type;
  /*const*/ std::array<FloatT, 3> position = {0.1, 0.2, 0.3};
  /*const*/ std::array<FloatT, 3> momentum = {10., 20., 30.};

  const auto pt = std::hypot(momentum[0], momentum[1]);
  const auto magField = 3.0;
  const auto charge = -1.0;

  helix.Initialize_VP(position.data(), momentum.data(), charge, magField);

  REQUIRE(helix.getRadius() == Catch::Approx(pt / (c_helix * magField)));
  REQUIRE(helix.getOmega() == Catch::Approx(charge / helix.getRadius()));
  REQUIRE(helix.getTanLambda() == Catch::Approx(momentum[2] / pt));
  // etc...
}

// TODO: add actually useful HelixClass tests, e.g.
// - Initialize with one of the methods and check whether the resulting Helix
// has the expected properties, like
//   - calculating different distances to the Helix and check whether the
//   conform to expecations
//   - Calculating momenta and see whether that gives the expected values
//   - essentially all non-trivial get functionality
// - Other usage examples as found in the "real-world"
// - Also add things that work slightly unexpectedly at the moment, e.g. things
// that are listed in https://github.com/iLCSoft/MarlinUtil/issues/24
