extern crate cgcamera;
extern crate cgmath;


use cgcamera::{CameraAttitude, Frustum, OrthographicCamera};
use cgmath::{Matrix4};


#[test]
fn test_unit_viewing_volume_projection_matrix_should_be_mirroring_matrix() {
    let left = -1.0;
    let right = 1.0;
    let bottom = -1.0;
    let top = 1.0;
    let near = 1.0;
    let far = -1.0;
    let frustum = Frustum::new(left, right, bottom, top, near, far);

    let origin = cgmath::vec3((0.0, 0.0, 5.0));
    let forward = cgmath::vec4((0.0, 0.0, -1.0, 0.0));
    let right = cgmath::vec4((1.0, 0.0, 0.0, 0.0));
    let up = cgmath::vec4((0.0, 1.0, 0.0, 0.0));
    let rotation_axis = cgmath::vec3((0.0, 0.0, -1.0));
    let attitude = CameraAttitude::new(origin, forward, right, up, rotation_axis);
    let camera = OrthographicCamera::new(frustum, attitude);

    let mirror_mat = Matrix4::new(
        1.0, 0.0,  0.0, 0.0,
        0.0, 1.0,  0.0, 0.0,
        0.0, 0.0, -1.0, 0.0,
        0.0, 0.0,  0.0, 1.0
    );

    assert_eq!(camera.proj_mat, mirror_mat);
}
