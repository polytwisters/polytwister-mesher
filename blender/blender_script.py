"""This script is meant to be run in Blender, not the standard Python interpreter.
"""
import argparse
import json
import math
import pathlib
import sys
import traceback
import warnings

from typing import Optional, NamedTuple

import bpy
import mathutils

EXPECTED_BLENDER_VERSION = (4, 1)

# Radius of polytwister's minimum containing sphere.
# 20cm seems like a reasonable diameter for a physical polytwister sculpture.
DEFAULT_SCALE = 20e-2 / 2


HDRI_PATH = pathlib.Path(__file__).resolve().parent.parent / "assets/studio_environment_2k.exr"


####################################################################################################
# Utils & wrappers


def deselect_all():
    bpy.ops.object.select_all(action="DESELECT")


def do_scale(amount):
    bpy.ops.object.origin_set(type="ORIGIN_CURSOR")
    bpy.ops.transform.resize(value=(amount, amount, amount))


def rotate_about_axis(axis, angle):
    # See https://stackoverflow.com/a/67697363.
    for area in bpy.context.screen.areas:
        if area.type == "VIEW_3D":
            view_3d = area
            break 
    else:
        raise RuntimeError("VIEW_3D area not found")

    bpy.ops.object.origin_set(type="ORIGIN_CURSOR")
    with bpy.context.temp_override(area=view_3d):
        bpy.ops.transform.rotate(
            value=angle,
            orient_axis=axis,
        )


def group_under_empty(parts):
    """Create an empty and group every object in parts as the child of
    that empty."""
    bpy.ops.object.empty_add(type="PLAIN_AXES")
    parent = bpy.context.object
    for part in parts:
        part.select_set(True)
    parent.select_set(True)
    bpy.ops.object.parent_set(type="OBJECT")
    deselect_all()
    bpy.context.view_layer.objects.active = parent
    parent.select_set(True)
    return parent


def rotation_to_point_to_origin(point):
    """Given a location of an object pointing along the X-axis, return a set of
    Euler angles that will rotate that object so it points at the origin."""
    # See https://blender.stackexchange.com/a/5220.
    direction = -mathutils.Vector(point)
    rotation_quaternion = direction.to_track_quat("-Z", "Y")
    return rotation_quaternion.to_euler()


def convert_spherical_to_cartesian(radius, latitude, longitude):
    """Convert radius, latitude, and longitude to Cartesian coordinates.
    Angles are in radians. The north pole/south pole axis is the Z axis,
    which looks vertical by default in Blender.

    Note the use of latitude (signed angle from equator) rather than
    zenith angle (angle from the north pole).
    """
    return (
        radius * math.cos(longitude) * math.cos(latitude),
        radius * math.sin(longitude) * math.cos(latitude),
        radius * math.sin(latitude),
    )


####################################################################################################
# Render setup


def set_up_camera(camera_longitude):
    """Add and return a camera object and set it to the primary camera of the scene. The input
    longitude is given in radians. A longitude of 0 is located on the positive X-axis, pi/2 on the
    positive Y-axis, etc.

    The camera is positioned about 15 degrees above the XY-plane and pointed directly at the origin.
    Its size in the viewport is also reduced.

    The camera is set to 85mm focal length, which is longer than Blender's default. This reduces
    the foreshortening while remaining realistic.
    """
    camera_distance = 1.0

    camera_latitude = math.radians(15)
    camera_location = convert_spherical_to_cartesian(
        camera_distance, camera_latitude, camera_longitude
    )
    bpy.ops.object.camera_add(
        location=camera_location,
        rotation=rotation_to_point_to_origin(camera_location)
    )
    camera = bpy.context.object
    camera.data.lens = 85
    camera.data.display_size = 0.1
    bpy.context.scene.camera = camera

    return camera


def set_up_environment(strength, image_path):
    """Set the world environment to an image, usually an HDRI environment. Its strength can be
    controlled. The environment does not appear in the result because of the transparent background.
    """
    world_node_tree = bpy.context.scene.world.node_tree
    nodes = world_node_tree.nodes
    nodes.clear()
    background_node = nodes.new(type="ShaderNodeBackground")
    background_node.inputs["Strength"].default_value = strength
    environment_node = nodes.new(type="ShaderNodeTexEnvironment")
    environment_node.image = bpy.data.images.load(image_path)
    output_node = nodes.new(type="ShaderNodeOutputWorld")
    links = world_node_tree.links
    links.new(environment_node.outputs["Color"], background_node.inputs["Color"])
    links.new(background_node.outputs["Background"], output_node.inputs["Surface"])


def set_up_lights(camera_longitude):
    """Add big, soft area lights. Hard shadows can give the impression of features that aren't
    really there, and aren't appropriate for presenting mathematical objects. However, we also must
    be careful to ensure nice contrast. If it is too washed out, it looks unattractive and inhibits
    depth perception.

    A bright key light highlights the object and a fill light ensures the shadows are not too dark.
    """
    # Latitude and longitude are given in degrees for readability.
    # Longitudes are relative to the camera: a longitude of 0 degrees
    # is directly behind the camera, 90 degrees is directly from the
    # right, -90 degrees is directly from the left. Latitudes are
    # absolute.
    light_specs = [
        # Key light illuminates most of the front of the object
        {"latitude": 10.0, "longitude": -80, "power": 300.0, "radius": 3.0},
        # Fill light gently illuminates the shadows left by the key light
        # Don't make this too strong, shadows are good
        {"latitude": 0.0, "longitude": 50.0, "power": 15.0, "radius": 3.0},
        # I used to have a back light here but it didn't work too well. The HDRI is sufficient for
        # preventing really dark areas.
    ]
    # If the lights are too bright, then this variable allows dimming all at once.
    power_multiplier = 1.0
    distance = 5.0

    for light_spec in light_specs:
        latitude = math.radians(light_spec["latitude"])
        longitude = camera_longitude + math.radians(light_spec["longitude"])
        location = convert_spherical_to_cartesian(
            distance, latitude, longitude
        )
        bpy.ops.object.light_add(
            type="AREA",
            radius=light_spec["radius"],
            location=location,
            rotation=rotation_to_point_to_origin(location)
        )
        bpy.context.object.data.energy = light_spec["power"] * power_multiplier


def set_transparent_background():
    """Ensure the environment texture does not show up in the video."""
    bpy.context.scene.render.film_transparent = True


def set_ambient_occlusion():
    """Enable "fast GI" (a type of ambient occlusion) that adds some darkness to the crevices of
    nonconvex polytwisters."""
    bpy.context.scene.cycles.use_fast_gi = True
    bpy.context.scene.world.light_settings.distance = 0.2


def set_image_size(resolution, camera):
    """Set the image size to resolution x resolution in pixels."""
    old_resolution = max([
        bpy.context.scene.render.resolution_x,
        bpy.context.scene.render.resolution_y,
    ])
    bpy.context.scene.render.resolution_x = resolution
    bpy.context.scene.render.resolution_y = resolution

    # Adjust the camera width to compensate for the change in resolution.
    # See https://blender.stackexchange.com/a/105805/154615.
    camera.data.sensor_width *= resolution / old_resolution


def set_render_engine():
    """Enable Cycles and hardware acceleration."""
    bpy.context.scene.render.engine = "CYCLES"
    bpy.context.scene.cycles.device = "GPU"


def set_sample_count(samples, preview_samples):
    """Adjust the Cycles sample count for the render and preview. The defaults are absurdly large
    (4096 for the render!)."""
    bpy.context.scene.cycles.samples = samples
    bpy.context.scene.cycles.preview_samples = preview_samples


def set_up_for_render(config):
    """Given a dictionary of render configuration settings, set the following:

    - Camera
    - Environment texture
    - Lighting
    - Ambient occlusion
    - Transparent background
    - Image resolution
    - Render engine
    - Sample count
    - Look
    """
    camera_longitude = math.radians(10)
    camera = set_up_camera(camera_longitude)

    set_up_environment(
        config.get("environment_strength", 0.3),
        config.get("environment_image", str(HDRI_PATH))
    )
    set_up_lights(camera_longitude)
    set_ambient_occlusion()
    set_transparent_background()
    set_image_size(config.get("resolution", 1080), camera)
    set_render_engine()
    set_sample_count(config.get("samples", 16), config.get("preview_samples", 4))


####################################################################################################
# Materials


class MaterialConfigs(NamedTuple):
    rings: dict
    strips: dict
    twisters_1: dict
    twisters_2: dict

    @classmethod
    def default(cls):
        ring_material_config = {
            "Base Color": [
                0.5,
                0.5,
                0.5,
                1
            ],
            "Roughness": 0.5,
        }

        strip_material_config = {
            "Base Color": [
                1.0,
                1.0,
                1.0,
                1
            ],
            "Roughness": 0.5,
        }

        twister_1_material_config = {
            "Base Color": [
                0.948,
                0.1,
                0.2,
                1
            ],
            "Roughness": 0.5,
        }
        twister_2_material_config = {
            "Base Color": [
                0.148,
                0.338,
                1.0,
                1
            ],
            "Roughness": 0.5,
        }
        return cls(
            rings=ring_material_config,
            strips=strip_material_config,
            twisters_1=twister_1_material_config,
            twisters_2=twister_2_material_config,
        )


def make_material_from_config(config):
    """Create a material on the active object and configure a Principled BSDF."""
    bpy.ops.material.new()
    material = bpy.data.materials[-1]

    if config is None:
        return material

    principled_bsdf = material.node_tree.nodes["Principled BSDF"]
    for key, value in config.items():
        if isinstance(value, list):
            value = tuple(value)
        principled_bsdf.inputs[key].default_value = value

    return material


materials_by_config = {}

def make_material_if_needed(config):
    """Memoized version of make_material_from_config."""
    # Stupid way of making the keys hashable.
    key = str(config)
    if key in materials_by_config:
        return materials_by_config[key]
    material = make_material_from_config(config)
    materials_by_config[key] = material
    return material


####################################################################################################
# PLY imports


def ply_is_empty(path):
    with open(path, mode='rb') as file:
        first_bytes = file.read(128)
    return first_bytes.split(b"\n")[2] == b"element vertex 0"


def import_ply(
    path: pathlib.Path,
    material_config: dict,
    frame_number: Optional[int] = None,
):
    # Blender creates an error importing an empty PLY file.
    if ply_is_empty(path):
        return

    deselect_all()
    bpy.ops.wm.ply_import(filepath=str(path))
    bpy.context.object.active_material = make_material_if_needed(material_config)

    bpy.ops.object.shade_smooth()
    bpy.ops.mesh.customdata_custom_splitnormals_clear()
    do_scale(DEFAULT_SCALE)

    if frame_number is not None:
        # To animate the sections, drivers are added so that the object appears only for its
        # assigned frame, both in the viewport and in the render.
        #
        # I decided not to use keyframes because the hide_viewport property can't be animated, which
        # is annoying. Drivers are fairly similar to keyframes in Blender and I haven't noticed any
        # performance issues in the viewport.
        driver = bpy.context.object.driver_add("hide_viewport").driver
        driver.type = "SCRIPTED"
        driver.expression = f"frame != {frame_number}"

        driver = bpy.context.object.driver_add("hide_render").driver
        driver.type = "SCRIPTED"
        driver.expression = f"frame != {frame_number}"

    return bpy.context.active_object


def import_plys(
    section_dir: pathlib.Path,
    material_configs: MaterialConfigs,
    frame_number: Optional[int] = None,
):
    """Import PLY files for a single section."""
    import_ply(
        section_dir / "rings.ply",
        frame_number=frame_number,
        material_config=material_configs.rings,
    )
    import_ply(
        section_dir / "strips.ply",
        frame_number=frame_number,
        material_config=material_configs.strips,
    )
    import_ply(
        section_dir / "twisters_1.ply",
        frame_number=frame_number,
        material_config=material_configs.twisters_1,
    )
    import_ply(
        section_dir / "twisters_2.ply",
        frame_number=frame_number,
        material_config=material_configs.twisters_2,
    )


def import_animation(root_dir: pathlib.Path, material_config: MaterialConfigs):
    section_dirs = []

    i = 0
    while True:
        section_dir = directory / f"section_{i:04}.ply"
        if not section_dir.exists():
            break
        section_dirs.append(section_dir)
        i += 1
    
    num_proper_frames = len(ring_paths)
    # One empty frame is added to the beginning and end of the animation. All frames in the middle
    # I call "proper frames."
    num_frames = num_proper_frames + 2
    bpy.context.scene.frame_end = num_frames

    for i in range(num_proper_frames):
        # +1 to convert 0-indexing to 1-indexing, another +1 for the initial empty frame.
        frame_number = i + 2
        import_plys(
            section_dirs[i],
            material_configs=material_config,
            frame_number=frame_number,
        )
    
    # To make things a bit more convenient when opening the .blend file interactively, navigate to
    # a frame where there is a visible mesh and align with the camera.
    bpy.context.scene.frame_set(num_proper_frames // 2)
    area = next(area for area in bpy.context.screen.areas if area.type == "VIEW_3D")
    area.spaces[0].region_3d.view_perspective = "CAMERA"


def is_animation_dir(root_dir: pathlib.Path) -> bool:
    return (root_dir / "section_0000").exists()


def main():
    major, minor, patch = bpy.app.version
    if (major, minor) != EXPECTED_BLENDER_VERSION:
        warnings.warn(
            "This script is tested with Blender "
            f"{EXPECTED_BLENDER_VERSION[0]}.{EXPECTED_BLENDER_VERSION[1]}, but you have "
            f"{major}.{minor}.{patch}. Errors may occur."
        )

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "dir",
        help="Input dir. If it contains a subdirectory named 'section_0000' it is imported as animation, otherwise a single section.",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="If provided, saves a .blend file to the given location.",
    )
    parser.add_argument(
        "-c",
        "--config",
        help="If provided, a JSON config string.",
    )

    argv = sys.argv
    for i, argument in enumerate(argv):
        if argument == "--":
            argv = argv[i + 1:]
            break
    else:
        argv = []
    args = parser.parse_args(argv)

    # Delete the default objects.
    bpy.ops.object.select_all(action="SELECT")
    bpy.ops.object.delete(use_global=False)

    directory = pathlib.Path(args.dir)

    if args.config is not None:
        config = json.loads(args.config)
    else:
        config = {}

    render_config = config.get("render", {})
    set_up_for_render(render_config)

    materials_config = MaterialConfigs.default()

    if is_animation_dir(directory):
        import_animation(directory, materials_config)
    else:
        import_plys(directory, materials_config)

    if args.output:
        # save_as_mainfile doesn't like relative paths, convert to absolute.
        path = pathlib.Path(args.output).resolve()
        # Pack resources
        bpy.ops.file.pack_all()
        # For info on relative_remap: https://blender.stackexchange.com/a/124861/154615
        bpy.ops.wm.save_as_mainfile(filepath=str(path), relative_remap=False)


if __name__ == "__main__":
    try:
        main()
    except Exception:
        traceback.print_exc()
        sys.exit(1)
