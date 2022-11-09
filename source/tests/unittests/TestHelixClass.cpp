#include "HelixClass.h"
#include "HelixClass_double.h"

#include "CLHEP/Units/PhysicalConstants.h"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_quantifiers.hpp>

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


TEMPLATE_LIST_TEST_CASE("getDistance-Initialize_VP",
                        "[helix-init][helix-distance]", HelixTypes) {
  TestType helix1, helix2;
  using FloatT = typename TestType::float_type;
  /*const*/ std::array<FloatT, 3> momentum1 = {10., 20., 30.};
  /*const*/ std::array<FloatT, 3> momentum2 = {70., 30., 10.};
  /*const*/ std::array<FloatT, 3> momentum;
  for (int i=0; i<3; i++) momentum[i] = momentum1[i] + momentum2[i];

  const auto magField = 3.0;
  const auto charge1 = -1.0;
  const auto charge2 = 1.0;

  FloatT vtx_momentum1[3];
  FloatT vtx_momentum2[3];
  FloatT vertex1[3];
  FloatT vertex2[3];

  SECTION("zero reference point") {
    // use ref. point in a "true" vertex for simplicity
    /*const*/ std::array<FloatT, 3> position = {0, 0, 0};

    helix1.Initialize_VP(position.data(), momentum1.data(), charge1, magField);
    helix2.Initialize_VP(position.data(), momentum2.data(), charge2, magField);

    FloatT dist1 = helix1.getDistanceToHelix(&helix2, vertex1, vtx_momentum1);
    FloatT dist2 = helix2.getDistanceToHelix(&helix1, vertex2, vtx_momentum2);

    FloatT eps = 0.1;

    // distance should be 0 as we initialized helices at the same point
    CHECK_THAT(0.0, Catch::Matchers::WithinAbs(dist1, eps));
    CHECK_THAT(0.0, Catch::Matchers::WithinAbs(dist2, eps));

    for (int i=0; i<3; i++) {
      CHECK_THAT(position[i], Catch::Matchers::WithinAbs(vertex1[i], eps));
      CHECK_THAT(position[i], Catch::Matchers::WithinAbs(vertex2[i], eps));

      CHECK_THAT(momentum[i], Catch::Matchers::WithinAbs(vtx_momentum1[i], eps));
      CHECK_THAT(momentum[i], Catch::Matchers::WithinAbs(vtx_momentum2[i], eps));
    }
  }

  SECTION("non-zero reference point") {
    // use ref. point in a "true" vertex for simplicity
    /*const*/ std::array<FloatT, 3> position = {214, 360, -700};

    helix1.Initialize_VP(position.data(), momentum1.data(), charge1, magField);
    helix2.Initialize_VP(position.data(), momentum2.data(), charge2, magField);

    FloatT dist1 = helix1.getDistanceToHelix(&helix2, vertex1, vtx_momentum1);
    FloatT dist2 = helix2.getDistanceToHelix(&helix1, vertex2, vtx_momentum2);

    FloatT eps = 0.1;

    CHECK_THAT(0.0, Catch::Matchers::WithinAbs(dist1, eps));
    CHECK_THAT(0.0, Catch::Matchers::WithinAbs(dist2, eps));

    for (int i=0; i<3; i++) {
      CHECK_THAT(position[i], Catch::Matchers::WithinAbs(vertex1[i], eps));
      CHECK_THAT(position[i], Catch::Matchers::WithinAbs(vertex2[i], eps));

      CHECK_THAT(momentum[i], Catch::Matchers::WithinAbs(vtx_momentum1[i], eps));
      CHECK_THAT(momentum[i], Catch::Matchers::WithinAbs(vtx_momentum2[i], eps));
    }
  }

}

TEMPLATE_LIST_TEST_CASE("getDistance-Initialize_Canonical",
                        "[helix-init][helix-distance]", HelixTypes) {
  TestType helix1, helix2;
  using FloatT = typename TestType::float_type;
  // "real-world" numbers (BSM event)
  /*const*/ std::array<FloatT, 3> momentum1 = {2.63e-01,-1.17e-01,-1.46e-01};
  /*const*/ std::array<FloatT, 3> momentum2 = {1.69e-01, 4.80e-02, 6.35e-01};
  /*const*/ std::array<FloatT, 3> momentum;
  for (int i=0; i<3; i++) momentum[i] = momentum1[i] + momentum2[i];

  const auto magField = 3.0;
  const auto charge1 = -1.0;
  const auto charge2 = 1.0;
  const auto pt1 = std::hypot(momentum1[0], momentum1[1]);
  const auto pt2 = std::hypot(momentum2[0], momentum2[1]);

  const auto omega1 = c_helix * magField * charge1 / pt1;
  const auto omega2 = c_helix * magField * charge2 / pt2;
  const auto tanLambda1 = momentum1[2] / pt1;
  const auto tanLambda2 = momentum2[2] / pt2;
  const auto phi01 = atan2( momentum1[1],momentum1[0] );
  const auto phi02 = atan2( momentum2[1],momentum2[0] );


  FloatT vtx_momentum1[3];
  FloatT vtx_momentum2[3];
  FloatT vertex1[3];
  FloatT vertex2[3];

  SECTION("no reference point given and helices passing through zero") {
    /*const*/ std::array<FloatT, 3> position = {0, 0, 0};

    // passing through ref. point
    const auto d0 = 1.0e-11;
    const auto z0 = 1.0e-11;

    helix1.Initialize_Canonical(phi01, d0, z0, omega1, tanLambda1, magField);
    helix2.Initialize_Canonical(phi02, -d0, z0, omega2, tanLambda2, magField);

    FloatT dist1 = helix1.getDistanceToHelix(&helix2, vertex1, vtx_momentum1);
    FloatT dist2 = helix2.getDistanceToHelix(&helix1, vertex2, vtx_momentum2);

    FloatT eps = 0.1;

    CHECK_THAT(0.0, Catch::Matchers::WithinAbs(dist1, eps));
    CHECK_THAT(0.0, Catch::Matchers::WithinAbs(dist2, eps));

    for (int i=0; i<3; i++) {
      CHECK_THAT(position[i], Catch::Matchers::WithinAbs(vertex1[i], eps));
      CHECK_THAT(position[i], Catch::Matchers::WithinAbs(vertex2[i], eps));

      CHECK_THAT(momentum[i], Catch::Matchers::WithinAbs(vtx_momentum1[i], eps));
      CHECK_THAT(momentum[i], Catch::Matchers::WithinAbs(vtx_momentum2[i], eps));
    }
  }

  SECTION("input zero reference point and zero d0, z0") {
    std::array<FloatT, 3> position = {0, 0, 0};

    // passing through ref. point (0,0,0)
    const auto d0 = 1.0e-11;
    const auto z0 = 1.0e-11;

    helix1.Initialize_Canonical(phi01, d0, z0, omega1, tanLambda1,
                                magField, position.data());
    helix2.Initialize_Canonical(phi02, -d0, z0, omega2, tanLambda2,
                                magField, position.data());

    FloatT dist1 = helix1.getDistanceToHelix(&helix2, vertex1, vtx_momentum1);
    FloatT dist2 = helix2.getDistanceToHelix(&helix1, vertex2, vtx_momentum2);

    FloatT eps = 0.1;

    CHECK_THAT(0.0, Catch::Matchers::WithinAbs(dist1, eps));
    CHECK_THAT(0.0, Catch::Matchers::WithinAbs(dist2, eps));

    for (int i=0; i<3; i++) {
      CHECK_THAT(position[i], Catch::Matchers::WithinAbs(vertex1[i], eps));
      CHECK_THAT(position[i], Catch::Matchers::WithinAbs(vertex2[i], eps));

      CHECK_THAT(momentum[i], Catch::Matchers::WithinAbs(vtx_momentum1[i], eps));
      CHECK_THAT(momentum[i], Catch::Matchers::WithinAbs(vtx_momentum2[i], eps));
    }
  }

  SECTION("input reference point and helices with small d0, z0") {
    // "real-world" numbers (BSM event)
    std::array<FloatT, 3> position = {2.42e+02,-3.58e+02, 1.72e+02}; // true vtx
    std::array<FloatT, 3> ref1 = {244.944, -359.603, 169.052};
    std::array<FloatT, 3> ref2 = {249.609, -356.381, 198.409};

    helix1.Initialize_Canonical(phi01, 0.103768967092, 1.39691245556, omega1,
                                tanLambda1, magField, ref1.data());
    helix2.Initialize_Canonical(phi02, 0.275440096855, -1.55072879791, omega2,
                                tanLambda2, magField, ref2.data());

    FloatT dist1 = helix1.getDistanceToHelix(&helix2, vertex1, vtx_momentum1);
    FloatT dist2 = helix2.getDistanceToHelix(&helix1, vertex2, vtx_momentum2);

    // distance should be small but nonzero
    CHECK_THAT(0.0, Catch::Matchers::WithinAbs(dist1, 2));
    CHECK_THAT(0.0, Catch::Matchers::WithinAbs(dist2, 2));

    for (int i=0; i<3; i++) {
      // 1 mm precision is good enough for very far vertices
      CHECK_THAT(position[i], Catch::Matchers::WithinAbs(vertex1[i], 1));
      CHECK_THAT(position[i], Catch::Matchers::WithinAbs(vertex2[i], 1));

      CHECK_THAT(momentum[i], Catch::Matchers::WithinAbs(vtx_momentum1[i], 0.1));
      CHECK_THAT(momentum[i], Catch::Matchers::WithinAbs(vtx_momentum2[i], 0.1));
    }
  }

  SECTION("input zero reference points and helices with large d0, z0") {
    // "real-world" numbers (BSM event)
    std::array<FloatT, 3> position = {2.42e+02,-3.58e+02, 1.72e+02};
    std::array<FloatT, 3> ref = {0, 0, 0};

    helix1.Initialize_Canonical(5.20e-01, -4.28e+02, 5.215e+01, 5.96e-03,
                                3.607e+00, magField, ref.data());
    helix2.Initialize_Canonical(1.27e+00, -9.89e+01, 3.709e+02, 3.64e-03,
                                5.050e-01, magField, ref.data());

    FloatT dist1 = helix1.getDistanceToHelix(&helix2, vertex1, vtx_momentum1);
    FloatT dist2 = helix2.getDistanceToHelix(&helix1, vertex2, vtx_momentum2);

    // distance should be small but nonzero
    CHECK_THAT(0.0, Catch::Matchers::WithinAbs(dist1, 2));
    CHECK_THAT(0.0, Catch::Matchers::WithinAbs(dist2, 2));

    for (int i=0; i<3; ++i) {
      // 1 mm precision is good enough for very far vertices
      CHECK_THAT(position[i], Catch::Matchers::WithinAbs(vertex1[i], 1));
      CHECK_THAT(position[i], Catch::Matchers::WithinAbs(vertex2[i], 1));

      CHECK_THAT(momentum[i], Catch::Matchers::WithinAbs(vtx_momentum1[i], 0.1));
      CHECK_THAT(momentum[i], Catch::Matchers::WithinAbs(vtx_momentum2[i], 0.1));
    }
  }
}
