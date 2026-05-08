# Simulation examples

## Simply supported I-beam three-point bending

`simply_supported_i_beam_bending.py` implements the problem from Figure 3:

- simple span `L = 4 m`, midpoint load `F = 10 kN`;
- steel with `E = 210 GPa`, `nu = 0.28`;
- I-section with `h = 0.2 m`, `bf = 0.1 m`, `tf = tw = 0.01 m`;
- self-weight ignored.

Run the local finite-element verification and regenerate the COMSOL Java Shell setup script:

```bash
python examples/simply_supported_i_beam_bending.py
```

Expected beam-element results are approximately:

- strong-axis second moment: `I = 2.292667e-05 m^4`;
- midpoint deflection: `w_mid = -2.769354e-03 m`;
- end rotations: `theta_left = -2.077016e-03 rad`, `theta_right = 2.077016e-03 rad`.

The generated `simply_supported_i_beam_bending.java` file contains COMSOL Java Shell commands for two GUI-checkable variants:

1. a beam-element model using a three-node line geometry;
2. a 3D solid-element model using an extruded I-section.

COMSOL entity numbering can differ between versions and model histories, so verify the generated boundary and point selections in the GUI before solving the two COMSOL studies.
