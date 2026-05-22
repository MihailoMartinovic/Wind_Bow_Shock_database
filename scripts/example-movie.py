import numpy as np
import pyvista as pv
def main():
    # Off-screen rendering is the most robust default for batch scripts.
    pl = pv.Plotter(off_screen=True) # documented Plotter(...,off_screen=True) 17
    pl.set_background("white")
    
    # Add something to see: a sphere.
    sphere = pv.Sphere(radius=1.0)
    pl.add_mesh(sphere, smooth_shading=True)

    # Prepare an animated object: a small marker sphere.
    marker = pv.Sphere(radius=0.15)
    marker_actor = pl.add_mesh(marker, smooth_shading=True)

    pl.camera_position = "yz"
    # quick view preset; you can also set a 3-tuple cam pos

    pl.open_movie("demo.mp4", framerate=30)

    n_frames = 120
    for i in range(n_frames):
        t = 2 * np.pi * i / n_frames
        # Move the marker around the origin
        marker_actor.SetPosition(2.5 * np.cos(t), 2.5 * np.sin(t), 0.0)
        pl.write_frame() # 24
    pl.close()
    
if __name__ == "__main__":
    main()
