import glob
import os
import xml.etree.ElementTree as ET

import numpy as np
import pyvista as pv
from numpy import ndarray, dtype
from tqdm import tqdm

# -------------------------------Settings----------------------------------------------
BASE_PATH = r"files/"
SINGLE_FILE_NAME = "0.xml"
# receptor_color = "#9400D3"
# receptor_color = "#000000"
receptor_color = "#8B0000"
ligand_color = "#FF8200"

pv.global_theme.render_points_as_spheres = True
pv.global_theme.allow_empty_mesh = True

class InteractionClass:
    def __init__(self, plotter):
        self.plotter = plotter

        self.pc_receptors = None
        self.pc_ligands = None
        self.points_ligands, self.points_receptors = create_ligands_receptor(BASE_PATH + SINGLE_FILE_NAME)
        self.line_actors = list()

        self.ligands_state = True
        self.receptors_state = True
        self.mito_only_state = False
        self.bounding_box_state = False

        self.add_ligands_receptors_to_plot()

        # Add key event to plotter
        self.plotter.add_key_event('c', self.activate_camera_tracker)
        self.plotter.add_key_event('l', self.hide_ligands)
        self.plotter.add_key_event('r', self.hide_receptors)
        self.plotter.add_key_event('o', self.get_camera_position)
        self.plotter.add_key_event('p', self.create_screenshot)
        self.plotter.add_key_event('t', self.reset_camera)
        self.plotter.add_key_event('h', self.print_help)
        self.plotter.add_key_event('m', self.show_only_mito)
        self.plotter.add_key_event('b', self.show_bounding_box)
        self.plotter.add_key_event('z', self.only_receptors)

    def activate_camera_tracker(self):
        self.plotter.add_camera_orientation_widget(animate=True, n_frames=20)

    def hide_ligands(self):
        if self.ligands_state:
            self.plotter.clear_actors()
            create_mitos(BASE_PATH + SINGLE_FILE_NAME, plotter=self.plotter)
            self.plotter.add_mesh(self.pc_receptors, render_points_as_spheres=True, color=receptor_color, point_size=3)
        else:
            self.plotter.clear_actors()
            create_mitos(BASE_PATH + SINGLE_FILE_NAME, plotter=self.plotter)
            self.add_ligands_receptors_to_plot()
        self.ligands_state = not self.ligands_state

    def hide_receptors(self):
        if self.receptors_state:
            self.plotter.clear_actors()
            create_mitos(BASE_PATH + SINGLE_FILE_NAME, plotter=self.plotter)
            self.plotter.add_mesh(self.pc_ligands, render_points_as_spheres=True, color=ligand_color, point_size=3)
        else:
            self.plotter.add_mesh(self.pc_receptors, render_points_as_spheres=True, color=receptor_color, point_size=3)
        self.receptors_state = not self.receptors_state

    def only_receptors(self):
        if self.receptors_state:
            self.plotter.clear_actors()
        else:
            self.plotter.add_mesh(self.pc_receptors, render_points_as_spheres=True, color=receptor_color, point_size=3)
        self.receptors_state = not self.receptors_state

    def show_only_mito(self):
        if self.mito_only_state:
            self.plotter.add_mesh(self.pc_receptors, render_points_as_spheres=True, color=receptor_color, point_size=3)
            self.plotter.add_mesh(self.pc_ligands, render_points_as_spheres=True, color=ligand_color, point_size=3)
        else:
            self.plotter.clear_actors()
            create_mitos(BASE_PATH + SINGLE_FILE_NAME, plotter=self.plotter)
        self.mito_only_state = not self.mito_only_state

    def get_camera_position(self):
        print(self.plotter.camera_position)
        print(self.plotter.camera.position)
        print(self.plotter.camera.azimuth)

    def create_screenshot(self):
        self.plotter.screenshot(
            filename=f"./output/screenshot_{os.path.basename(os.path.normpath(BASE_PATH))}_{SINGLE_FILE_NAME.replace('.xml', '')}.png",
            transparent_background=True)

    def reset_camera(self):
        self.plotter.camera_position = [(43527, 97930, 5000),
                                        (43527, 2250, 5000),
                                        (0.0, 0.0, 1.0)]
        self.plotter.reset_camera_clipping_range()

    def add_ligands_receptors_to_plot(self):
        self.pc_receptors = pv.PolyData(self.points_receptors)
        self.pc_ligands = pv.PolyData(self.points_ligands)

        self.plotter.add_mesh(self.pc_receptors, render_points_as_spheres=True, color=receptor_color, point_size=3)
        self.plotter.add_mesh(self.pc_ligands, render_points_as_spheres=True, color=ligand_color, point_size=3)

    @staticmethod
    def print_help():
        print("H -- Help\nL -- Hide only Ligands\nR -- Hide only Receptors\nM -- Only Mitos\nT -- Reset Camera (press "
              "arrow key afterwards)\nP -- Screenshot (output/screenshot_filename\nO -- Print current camera position\n"
              "B -- Toggle Bounding Box\nZ -- Show only Receptors")

    def show_bounding_box(self):
        if not self.bounding_box_state:
            cell_x_max = 87053.77778
            cell_y_max = 4500
            cell_z_max = 10000
            
            line1 = pv.Line((0, 0, 0), (0, cell_y_max, 0))
            line2 = pv.Line((0, 0, 0), (0, 0, cell_z_max))
            line3 = pv.Line((0, 0, 0), (cell_x_max, 0, 0))

            line4 = pv.Line((0, cell_y_max, 0), (cell_x_max, cell_y_max, 0))
            line5 = pv.Line((0, cell_y_max, 0), (0, cell_y_max, cell_z_max))

            line6 = pv.Line((0, 0, cell_z_max), (cell_x_max, 0, cell_z_max))
            line7 = pv.Line((0, 0, cell_z_max), (0, cell_y_max, cell_z_max))

            line8 = pv.Line((0, cell_y_max, cell_z_max), (cell_x_max, cell_y_max, cell_z_max))

            line9 = pv.Line((cell_x_max, cell_y_max, cell_z_max), (cell_x_max, cell_y_max, 0))
            line10 = pv.Line((cell_x_max, cell_y_max, cell_z_max), (cell_x_max, 0, cell_z_max))

            line11 = pv.Line((cell_x_max, 0, 0), (cell_x_max, cell_y_max, 0))
            line12 = pv.Line((cell_x_max, 0, 0), (cell_x_max, 0, cell_z_max))

            lines = [line1, line2, line3, line4, line5, line6, line7, line8, line9, line10, line11, line12]
            for line in lines:
                self.line_actors.append(self.plotter.add_mesh(line, color="black", line_width=5))
        else:
            for actor in self.line_actors:
                self.plotter.remove_actor(actor)
        self.bounding_box_state = not self.bounding_box_state


def create_mitos(path: str, plotter: pv.Plotter) -> None:
    tree = ET.parse(path)
    root = tree.getroot()

    for xagent in tqdm(root.findall('xagent'), desc="Creating Mitos"):
        if xagent.find('name').text == "Mito":
            xmin = float(xagent.find('x1').text)
            ymin = float(xagent.find('y1').text)
            zmin = float(xagent.find('z1').text)

            xmax = float(xagent.find('x8').text)
            ymax = float(xagent.find('y8').text)
            zmax = float(xagent.find('z8').text)

            box = pv.Box(bounds=[zmin, zmax, ymin, ymax, xmin, xmax])

            plotter.add_mesh(box, color="lightblue", show_edges=True)

        else:
            continue


def create_ligands_receptor(path: str) -> tuple[ndarray[float, dtype[float]], ndarray[float, dtype[float]]]:
    tree = ET.parse(path)
    root = tree.getroot()

    points_ligands = list()
    points_receptors = list()
    for xagent in tqdm(root.findall('xagent'), desc="Creating ligands and receptors"):
        if xagent.find('name').text == "Ligand":
            x = float(xagent.find('x').text)
            y = float(xagent.find('y').text)
            z = float(xagent.find('z').text)

            # Transform Coordinates to make displaying them easier
            points_ligands.append([z, y, x])

        if xagent.find('name').text == "Receptor":
            x = float(xagent.find('x').text)
            y = float(xagent.find('y').text) - 0.1
            z = float(xagent.find('z').text)

            # Transform Coordinates to make displaying them easier
            points_receptors.append([z, y, x])

    points_ligands = np.asarray(points_ligands, dtype=np.float32)
    points_receptors = np.asarray(points_receptors, dtype=np.float32)

    return points_ligands, points_receptors


if __name__ == "__main__":
    plotter = pv.Plotter(image_scale=2, window_size=[1900, 1000])
    plotter.camera_position = [(43527, 97930, 5000),
                               (43527, 2250, 5000),
                               (0.0, 0.0, 1.0)]
    tracker = InteractionClass(plotter)

    create_mitos(BASE_PATH + SINGLE_FILE_NAME, plotter=plotter)

    plotter.reset_camera_clipping_range()
    plotter.show()
