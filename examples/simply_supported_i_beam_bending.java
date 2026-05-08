// Simply supported I-beam bending simulation: beam and solid FE variants.
// Paste into COMSOL Java Shell in a new 3D model, then verify selections.
model.modelPath(System.getProperty("user.dir"));
model.label("simply_supported_i_beam_three_point_bending.mph");

model.param().set("L", "4.0[m]", "Beam span");
model.param().set("F0", "10000.0[N]", "Downward center point load");
model.param().set("E0", "210000000000.0[Pa]", "Steel Young's modulus");
model.param().set("nu0", "0.28", "Steel Poisson ratio");
model.param().set("h", "0.2[m]", "I-section height");
model.param().set("bf", "0.1[m]", "Flange width");
model.param().set("tf", "0.01[m]", "Flange thickness");
model.param().set("tw", "0.01[m]", "Web thickness");
model.param().set("Asec", "0.0038[m^2]", "I-section area");
model.param().set("Iyy", "2.29266666667e-05[m^4]", "Strong-axis second moment of area");
model.param().set("G0", "82031250000[Pa]", "Steel shear modulus");

// ------------------------------------------------------------------
// Model 1: Beam elements on a 4 m line, with an explicit midpoint.
// ------------------------------------------------------------------
model.component().create("beamComp", true);
model.component("beamComp").geom().create("geomB", 3);
model.component("beamComp").geom("geomB").lengthUnit("m");
model.component("beamComp").geom("geomB").create("pt1", "Point");
model.component("beamComp").geom("geomB").feature("pt1").set("p", new String[]{"0", "0", "0"});
model.component("beamComp").geom("geomB").create("pt2", "Point");
model.component("beamComp").geom("geomB").feature("pt2").set("p", new String[]{"L/2", "0", "0"});
model.component("beamComp").geom("geomB").create("pt3", "Point");
model.component("beamComp").geom("geomB").feature("pt3").set("p", new String[]{"L", "0", "0"});
model.component("beamComp").geom("geomB").create("e1", "LineSegment");
model.component("beamComp").geom("geomB").feature("e1").selection("point1").set("pt1", 1);
model.component("beamComp").geom("geomB").feature("e1").selection("point2").set("pt2", 1);
model.component("beamComp").geom("geomB").create("e2", "LineSegment");
model.component("beamComp").geom("geomB").feature("e2").selection("point1").set("pt2", 1);
model.component("beamComp").geom("geomB").feature("e2").selection("point2").set("pt3", 1);
model.component("beamComp").geom("geomB").run();
model.component("beamComp").material().create("matSteelB", "Common");
model.component("beamComp").material("matSteelB").propertyGroup("def").set("youngsmodulus", "E0");
model.component("beamComp").material("matSteelB").propertyGroup("def").set("poissonsratio", "nu0");
model.component("beamComp").physics().create("beam", "Beam", "geomB");
// In the Beam interface, set cross-section data to: A=Asec, Iyy=Iyy, Izz as required by local axes.
// Apply: pinned support at x=0, roller support at x=L, point load Fy=-F0 at x=L/2.
model.component("beamComp").mesh().create("meshB");
model.component("beamComp").mesh("meshB").autoMeshSize(1);
model.component("beamComp").mesh("meshB").run();
model.study().create("stdBeam");
model.study("stdBeam").create("stat", "Stationary");

// ------------------------------------------------------------------
// Model 2: 3D solid elements. Build the I-section in the yz-plane and
// extrude it along x. Verify end-face and load-face selections in GUI.
// ------------------------------------------------------------------
model.component().create("solidComp", true);
model.component("solidComp").geom().create("geomS", 3);
model.component("solidComp").geom("geomS").lengthUnit("m");
model.component("solidComp").geom("geomS").create("wp1", "WorkPlane");
model.component("solidComp").geom("geomS").feature("wp1").set("quickplane", "yz");
model.component("solidComp").geom("geomS").feature("wp1").geom().create("topFlange", "Rectangle");
model.component("solidComp").geom("geomS").feature("wp1").geom().feature("topFlange").set("size", new String[]{"bf", "tf"});
model.component("solidComp").geom("geomS").feature("wp1").geom().feature("topFlange").set("pos", new String[]{"-bf/2", "h/2-tf"});
model.component("solidComp").geom("geomS").feature("wp1").geom().create("web", "Rectangle");
model.component("solidComp").geom("geomS").feature("wp1").geom().feature("web").set("size", new String[]{"tw", "h-2*tf"});
model.component("solidComp").geom("geomS").feature("wp1").geom().feature("web").set("pos", new String[]{"-tw/2", "-h/2+tf"});
model.component("solidComp").geom("geomS").feature("wp1").geom().create("botFlange", "Rectangle");
model.component("solidComp").geom("geomS").feature("wp1").geom().feature("botFlange").set("size", new String[]{"bf", "tf"});
model.component("solidComp").geom("geomS").feature("wp1").geom().feature("botFlange").set("pos", new String[]{"-bf/2", "-h/2"});
model.component("solidComp").geom("geomS").create("ext1", "Extrude");
model.component("solidComp").geom("geomS").feature("ext1").selection("input").set("wp1");
model.component("solidComp").geom("geomS").feature("ext1").set("distance", "L");
model.component("solidComp").geom("geomS").run();
model.component("solidComp").material().create("matSteelS", "Common");
model.component("solidComp").material("matSteelS").propertyGroup("def").set("youngsmodulus", "E0");
model.component("solidComp").material("matSteelS").propertyGroup("def").set("poissonsratio", "nu0");
model.component("solidComp").physics().create("solid", "SolidMechanics", "geomS");
// Apply simple supports to the two end faces and a distributed load over a small
// top-flange patch centered at x=L/2 with total vertical force -F0.
model.component("solidComp").mesh().create("meshS");
model.component("solidComp").mesh("meshS").autoMeshSize(3);
model.component("solidComp").mesh("meshS").run();
model.study().create("stdSolid");
model.study("stdSolid").create("stat", "Stationary");

// Suggested result checks:
// Beam midspan vertical displacement should be about -2.769e-3 m.
// Support rotations should be about +/-2.077e-3 rad.
