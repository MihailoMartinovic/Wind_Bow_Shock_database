import numpy as np
import pyvista as pv


def revolve_r_theta_to_surface(theta, r_theta, n_phi=72):
    theta = np.asarray(theta, float)
    r_theta = np.asarray(r_theta, float)

    phi = np.linspace(0.0, 2.0 * np.pi, n_phi, endpoint=False)
    th, ph = np.meshgrid(theta, phi, indexing="ij")

    r = r_theta[:, None] * np.ones_like(ph)

    x = r * np.cos(th)
    rho = r * np.sin(th)
    y = rho * np.cos(ph)
    z = rho * np.sin(ph)

    grid = pv.StructuredGrid(x, y, z)
    surf = grid.extract_surface(algorithm="dataset_surface").triangulate()
    return surf


def polyline_tube(points_xyz, radius=0.05, n_sides=16):
    poly = pv.lines_from_points(points_xyz)
    return poly.tube(radius=radius, n_sides=n_sides)


def main():
    pl = pv.Plotter(off_screen=True, window_size=[1600, 1200])
    pl.set_background("white")

    earth = pv.Sphere(radius=1.0)
    pl.add_mesh(earth, color="lightblue", smooth_shading=True)

    # Back off from the theta -> pi singular tail
    theta = np.linspace(0.0, 0.95 * np.pi, 300)

    r0_mp = 10.0
    alpha_mp = 0.7
    r_mp = r0_mp * (2.0 / (1.0 + np.cos(theta))) ** alpha_mp
    r_mp = np.clip(r_mp, None, 80.0)

    r0_bs = 15.0
    e = 0.9
    r_bs = r0_bs * (1.0 + e) / (1.0 + e * np.cos(theta))
    r_bs = np.clip(r_bs, None, 120.0)

    print("r_mp min/max:", np.min(r_mp), np.max(r_mp))
    print("r_bs min/max:", np.min(r_bs), np.max(r_bs))

    mp_surf = revolve_r_theta_to_surface(theta, r_mp, n_phi=96)
    bs_surf = revolve_r_theta_to_surface(theta, r_bs, n_phi=96)

    print("mp bounds:", mp_surf.bounds)
    print("bs bounds:", bs_surf.bounds)

    pl.add_mesh(mp_surf, color="dodgerblue", opacity=0.25)
    pl.add_mesh(bs_surf, color="orange", opacity=0.20)

    t = np.linspace(0, 10 * np.pi, 800)
    field_pts = np.c_[5 * np.cos(t), 2.5 * np.sin(t), 0.2 * t - 3]
    field_tube = polyline_tube(field_pts, radius=0.07)
    pl.add_mesh(field_tube, color="purple")

    tt = np.linspace(0, 2 * np.pi, 400)
    traj_pts = np.c_[20 * np.cos(tt) - 20, 8 * np.sin(tt), 4 * np.sin(2 * tt)]
    traj_tube = polyline_tube(traj_pts, radius=0.12)
    pl.add_mesh(traj_tube, color="black")

    sc_marker = pv.Sphere(radius=0.5)
    sc_actor = pl.add_mesh(sc_marker, color="red", smooth_shading=True)

    pl.camera_position = [
        (50, 30, 25),
        (-10, 0, 0),
        (0, 0, 1),
    ]

    pl.show(auto_close=False)
    pl.screenshot("debug_scene.png")

    pl.open_movie("magnetosphere_proto.mp4", framerate=30)

    n_frames = 240
    for i in range(n_frames):
        k = int((i / (n_frames - 1)) * (len(traj_pts) - 1))
        sc_actor.SetPosition(*traj_pts[k])
        pl.write_frame()

    pl.close()


if __name__ == "__main__":
    main()
