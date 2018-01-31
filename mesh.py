class MeshGrid

    # Specify boundary cells
    cells_top
    cells_bottom
    cells_left
    cells_right

    # Specify interior cells
    cells_interior

    construct_grid_from_geometry()
    visualize()
    visualize_colormap(u or T or p or etc..)

class MeshCell

    x, y, ID

    rho, u, T, p, etc...
    rho_old, u_old, T_old, p_old, etc...

    get_left_neighbor()
    get_right_neighbor()
    get_top_neighbor()
    get_bottom_neighbor()