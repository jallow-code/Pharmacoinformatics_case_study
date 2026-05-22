# PyMOL Presentation Rendering

Use the PyMOL scene setup files below to open each figure with the project styling already applied:

- `scripts/presentation_mutation_scenes.pml`
- `scripts/presentation_wt_panel.pml`
- `scripts/presentation_c326y_panel.pml`
- `scripts/presentation_y64n_panel.pml`
- `scripts/presentation_validation_pose_overlay.pml`
- `scripts/presentation_validation_box_overview.pml`
- `scripts/presentation_c326y_top_hits_overlay.pml`

## Open a scene

From the project root:

```bash
/home/jallow/miniforge3/bin/pymol scripts/presentation_wt_panel.pml
```

Replace the file name with the scene you want.

For the three mutation panels, use the combined mutation scene file:

```bash
/home/jallow/miniforge3/bin/pymol scripts/presentation_mutation_scenes.pml
```

Then switch between the stored scenes with:

```pml
scene WT, recall
scene C326Y, recall
scene Y64N, recall
```

## Adjust the angle

Inside PyMOL:

- left drag: rotate
- middle drag or wheel: zoom
- right drag: translate

The scripts already apply:

- white background
- orthoscopic projection
- shadows on
- ambient occlusion on
- presentation colors
- consistent cartoon and stick styling

## Render the final image

After you are happy with the view, run in the PyMOL command line:

```pml
ray 2400, 1800
png report/figures/presentation/my_figure.png, dpi=300
```

For the validation pose overlay and the C326Y top-hits overlay, use:

```pml
ray 3200, 2400
png report/figures/presentation/c326y_top_hits_overlay.png, dpi=300
```

## Keep the three mutation panels comparable

For the `WT`, `C326Y`, and `Y64N` panels:

1. open `presentation_mutation_scenes.pml`
2. recall `WT`
3. choose the angle
4. run `get_view`
5. render WT
6. recall `C326Y` and paste the same view with `set_view`
7. render C326Y
8. recall `Y64N` and apply the same `set_view`
9. render Y64N

Because the structures are pre-aligned, this is the cleanest way to keep the 1 row, 3 column panel matched.
