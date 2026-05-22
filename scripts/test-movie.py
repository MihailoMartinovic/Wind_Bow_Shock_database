import pyvista as pv

pl = pv.Plotter(off_screen=True, window_size=[800, 600])
pl.add_mesh(pv.Cube(), color="red")
pl.show(auto_close=False)          # initialize the scene
pl.screenshot("cube.png")
print("bounds:", pl.bounds)
print("camera:", pl.camera_position)
pl.close()
