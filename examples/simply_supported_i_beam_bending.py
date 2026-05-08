"""Finite-element simulation for a simply supported I-beam under a center load.

Problem statement represented by this script:
- Span: 4 m, simple supports at both ends.
- Center point load: 10 kN downward at x = 2 m.
- Steel: E = 210 GPa, Poisson ratio = 0.28, self-weight ignored.
- I-section: total height 0.2 m, flange width 0.1 m, flange/web thickness 0.01 m.

The script contains two complementary models:
1. A runnable Euler-Bernoulli beam finite-element model using two cubic beam
   elements. For this load case it reproduces the textbook three-point-bending
   solution exactly at the nodes.
2. A COMSOL Java Shell script generator that creates beam-element and solid-
   element model setup commands for the same geometry and boundary conditions.

Run locally with:
    python examples/simply_supported_i_beam_bending.py
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import textwrap


@dataclass(frozen=True)
class BeamInput:
    span_m: float = 4.0
    load_n: float = 10_000.0
    young_pa: float = 210e9
    poisson: float = 0.28
    section_height_m: float = 0.2
    flange_width_m: float = 0.1
    flange_thickness_m: float = 0.01
    web_thickness_m: float = 0.01

    @property
    def web_height_m(self) -> float:
        return self.section_height_m - 2.0 * self.flange_thickness_m

    @property
    def area_m2(self) -> float:
        return (
            2.0 * self.flange_width_m * self.flange_thickness_m
            + self.web_thickness_m * self.web_height_m
        )

    @property
    def second_moment_m4(self) -> float:
        """Strong-axis second moment of area of the I-section."""
        h = self.section_height_m
        b = self.flange_width_m
        tf = self.flange_thickness_m
        tw = self.web_thickness_m
        hw = self.web_height_m
        flange_centroid_offset = h / 2.0 - tf / 2.0
        flange_i = b * tf**3 / 12.0 + b * tf * flange_centroid_offset**2
        web_i = tw * hw**3 / 12.0
        return 2.0 * flange_i + web_i

    @property
    def shear_modulus_pa(self) -> float:
        return self.young_pa / (2.0 * (1.0 + self.poisson))


def euler_bernoulli_element_stiffness(ei: float, length: float) -> list[list[float]]:
    """Local bending stiffness matrix for DOFs [w1, theta1, w2, theta2]."""
    le = length
    scale = ei / le**3
    base = [
        [12.0, 6.0 * le, -12.0, 6.0 * le],
        [6.0 * le, 4.0 * le**2, -6.0 * le, 2.0 * le**2],
        [-12.0, -6.0 * le, 12.0, -6.0 * le],
        [6.0 * le, 2.0 * le**2, -6.0 * le, 4.0 * le**2],
    ]
    return [[scale * value for value in row] for row in base]


def solve_linear_system(matrix: list[list[float]], vector: list[float]) -> list[float]:
    """Solve Ax=b by Gaussian elimination with partial pivoting."""
    size = len(vector)
    augmented = [row[:] + [rhs] for row, rhs in zip(matrix, vector, strict=True)]

    for pivot_index in range(size):
        pivot_row = max(
            range(pivot_index, size), key=lambda row: abs(augmented[row][pivot_index])
        )
        if abs(augmented[pivot_row][pivot_index]) < 1e-30:
            raise ValueError("Singular stiffness matrix after applying constraints.")
        if pivot_row != pivot_index:
            augmented[pivot_index], augmented[pivot_row] = (
                augmented[pivot_row],
                augmented[pivot_index],
            )

        pivot = augmented[pivot_index][pivot_index]
        for col in range(pivot_index, size + 1):
            augmented[pivot_index][col] /= pivot

        for row in range(size):
            if row == pivot_index:
                continue
            factor = augmented[row][pivot_index]
            if factor == 0.0:
                continue
            for col in range(pivot_index, size + 1):
                augmented[row][col] -= factor * augmented[pivot_index][col]

    return [augmented[row][size] for row in range(size)]


def solve_two_element_beam(
    problem: BeamInput,
) -> dict[str, float | list[float] | list[list[float]]]:
    """Solve the three-node beam-element model.

    Nodes are placed at x = 0, L/2, L so that the point load is applied exactly
    at a finite-element node. Vertical translations are constrained at the two
    supports; all rotations remain free.
    """
    node_count = 3
    dof_per_node = 2
    total_dofs = node_count * dof_per_node
    element_length = problem.span_m / 2.0
    ei = problem.young_pa * problem.second_moment_m4

    stiffness = [[0.0 for _ in range(total_dofs)] for _ in range(total_dofs)]
    for element_index in range(2):
        element_dofs = [
            2 * element_index,
            2 * element_index + 1,
            2 * (element_index + 1),
            2 * (element_index + 1) + 1,
        ]
        ke = euler_bernoulli_element_stiffness(ei, element_length)
        for local_row, global_row in enumerate(element_dofs):
            for local_col, global_col in enumerate(element_dofs):
                stiffness[global_row][global_col] += ke[local_row][local_col]

    force = [0.0 for _ in range(total_dofs)]
    force[2] = -problem.load_n

    constrained = {0, 4}  # w(left)=0, w(right)=0
    free = [dof for dof in range(total_dofs) if dof not in constrained]

    reduced_stiffness = [[stiffness[row][col] for col in free] for row in free]
    reduced_force = [force[row] for row in free]
    reduced_displacement = solve_linear_system(reduced_stiffness, reduced_force)

    displacement = [0.0 for _ in range(total_dofs)]
    for dof, value in zip(free, reduced_displacement, strict=True):
        displacement[dof] = value

    reactions = [
        sum(stiffness[row][col] * displacement[col] for col in range(total_dofs))
        - force[row]
        for row in range(total_dofs)
    ]

    analytical_midspan_deflection = -problem.load_n * problem.span_m**3 / (
        48.0 * problem.young_pa * problem.second_moment_m4
    )
    analytical_end_rotation = problem.load_n * problem.span_m**2 / (
        16.0 * problem.young_pa * problem.second_moment_m4
    )

    return {
        "stiffness": stiffness,
        "force": force,
        "displacement": displacement,
        "reactions": reactions,
        "midspan_deflection_m": displacement[2],
        "left_rotation_rad": displacement[1],
        "right_rotation_rad": displacement[5],
        "analytical_midspan_deflection_m": analytical_midspan_deflection,
        "analytical_end_rotation_rad": analytical_end_rotation,
    }


def build_comsol_java_shell_script(problem: BeamInput) -> str:
    """Return COMSOL Java Shell commands for beam and solid FE models.

    The generated commands are intentionally parameterized and GUI-friendly: paste
    them into an empty COMSOL model's Java Shell, then inspect selections and run
    the studies. Boundary/point selections can vary between COMSOL versions, so
    comments mark the places that should be verified in the GUI before solving.
    """
    return textwrap.dedent(
        f"""
        // Simply supported I-beam bending simulation: beam and solid FE variants.
        // Paste into COMSOL Java Shell in a new 3D model, then verify selections.
        model.modelPath(System.getProperty("user.dir"));
        model.label("simply_supported_i_beam_three_point_bending.mph");

        model.param().set("L", "{problem.span_m}[m]", "Beam span");
        model.param().set("F0", "{problem.load_n}[N]", "Downward center point load");
        model.param().set("E0", "{problem.young_pa}[Pa]", "Steel Young's modulus");
        model.param().set("nu0", "{problem.poisson}", "Steel Poisson ratio");
        model.param().set("h", "{problem.section_height_m}[m]", "I-section height");
        model.param().set("bf", "{problem.flange_width_m}[m]", "Flange width");
        model.param().set("tf", "{problem.flange_thickness_m}[m]", "Flange thickness");
        model.param().set("tw", "{problem.web_thickness_m}[m]", "Web thickness");
        model.param().set("Asec", "{problem.area_m2:.12g}[m^2]", "I-section area");
        model.param().set("Iyy", "{problem.second_moment_m4:.12g}[m^4]", "Strong-axis second moment of area");
        model.param().set("G0", "{problem.shear_modulus_pa:.12g}[Pa]", "Steel shear modulus");

        // ------------------------------------------------------------------
        // Model 1: Beam elements on a 4 m line, with an explicit midpoint.
        // ------------------------------------------------------------------
        model.component().create("beamComp", true);
        model.component("beamComp").geom().create("geomB", 3);
        model.component("beamComp").geom("geomB").lengthUnit("m");
        model.component("beamComp").geom("geomB").create("pt1", "Point");
        model.component("beamComp").geom("geomB").feature("pt1").set("p", new String[]{{"0", "0", "0"}});
        model.component("beamComp").geom("geomB").create("pt2", "Point");
        model.component("beamComp").geom("geomB").feature("pt2").set("p", new String[]{{"L/2", "0", "0"}});
        model.component("beamComp").geom("geomB").create("pt3", "Point");
        model.component("beamComp").geom("geomB").feature("pt3").set("p", new String[]{{"L", "0", "0"}});
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
        model.component("solidComp").geom("geomS").feature("wp1").geom().feature("topFlange").set("size", new String[]{{"bf", "tf"}});
        model.component("solidComp").geom("geomS").feature("wp1").geom().feature("topFlange").set("pos", new String[]{{"-bf/2", "h/2-tf"}});
        model.component("solidComp").geom("geomS").feature("wp1").geom().create("web", "Rectangle");
        model.component("solidComp").geom("geomS").feature("wp1").geom().feature("web").set("size", new String[]{{"tw", "h-2*tf"}});
        model.component("solidComp").geom("geomS").feature("wp1").geom().feature("web").set("pos", new String[]{{"-tw/2", "-h/2+tf"}});
        model.component("solidComp").geom("geomS").feature("wp1").geom().create("botFlange", "Rectangle");
        model.component("solidComp").geom("geomS").feature("wp1").geom().feature("botFlange").set("size", new String[]{{"bf", "tf"}});
        model.component("solidComp").geom("geomS").feature("wp1").geom().feature("botFlange").set("pos", new String[]{{"-bf/2", "-h/2"}});
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
        """
    ).strip() + "\n"


def write_comsol_script(path: Path, problem: BeamInput) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(build_comsol_java_shell_script(problem), encoding="utf-8")


def main() -> None:
    problem = BeamInput()
    result = solve_two_element_beam(problem)
    output_script = Path(__file__).with_suffix(".java")
    write_comsol_script(output_script, problem)

    print("Simply supported I-beam three-point-bending simulation")
    print(f"A = {problem.area_m2:.6e} m^2")
    print(f"I = {problem.second_moment_m4:.6e} m^4")
    print(f"G = {problem.shear_modulus_pa:.6e} Pa")
    print(f"w_mid(FE) = {result['midspan_deflection_m']:.6e} m")
    print(f"theta_left(FE) = {result['left_rotation_rad']:.6e} rad")
    print(f"theta_right(FE) = {result['right_rotation_rad']:.6e} rad")
    print(f"w_mid(analytical) = {result['analytical_midspan_deflection_m']:.6e} m")
    print(f"theta_edge(analytical) = {result['analytical_end_rotation_rad']:.6e} rad")
    print(f"COMSOL Java Shell setup script written to: {output_script}")


if __name__ == "__main__":
    main()
