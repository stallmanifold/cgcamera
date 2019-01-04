use cgcamera::{CameraAttitude, Frustum, OrthographicCamera};
use cgmath::{Vector3, Matrix4};


fn orthographic_mirror_matrix_camera_model() -> OrthographicCamera {
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

    camera
}

fn mirror_mat() -> Matrix4 {
    Matrix4::new(
        1.0, 0.0,  0.0, 0.0,
        0.0, 1.0,  0.0, 0.0,
        0.0, 0.0, -1.0, 0.0,
        0.0, 0.0,  0.0, 1.0
    )
}

#[test]
fn test_unit_viewing_volume_projection_matrix_should_be_mirroring_matrix() {
    let camera = orthographic_mirror_matrix_camera_model();
    let proj_mat = mirror_mat();

    assert_eq!(camera.proj_mat, proj_mat);
}

#[test]
fn test_orthographic_camera_mirror_mat_trans_mat() {
    let camera = orthographic_mirror_matrix_camera_model();
    let trans_mat = Matrix4::new(
        1.0, 0.0,  0.0, 0.0,
        0.0, 1.0,  0.0, 0.0,
        0.0, 0.0,  1.0, 0.0,
        0.0, 0.0, -5.0, 1.0
    );

    assert_eq!(camera.trans_mat, trans_mat);
}

#[test]
fn test_orthographic_camera_mirror_mat_rot_mat() {
    let camera = orthographic_mirror_matrix_camera_model();
    let rot_mat = Matrix4::new(
        1.0, 0.0, 0.0, 0.0,
        0.0, 1.0, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, 1.0
    );

    assert_eq!(camera.rot_mat, rot_mat);
}

#[test]
fn test_orthographic_camera_mirror_mat_view_mat() {
    let camera = orthographic_mirror_matrix_camera_model();
    let view_mat = Matrix4::new(
        1.0, 0.0,  0.0,  0.0,
        0.0, 1.0,  0.0,  0.0,
        0.0, 0.0,  1.0,  0.0,
        0.0, 0.0, -5.0,  1.0
    );

    assert_eq!(camera.view_mat, view_mat);
}

fn orthographic_camera_z_axis() -> OrthographicCamera {
    let left = -4.0;
    let right = 4.0;
    let bottom = -4.0;
    let top = 4.0;
    let near = 4.0;
    let far = -4.0;
    let frustum = Frustum::new(left, right, bottom, top, near, far);

    let origin = cgmath::vec3((0.0, 0.0, 5.0));
    let forward = cgmath::vec4((0.0, 0.0, 1.0, 0.0));
    let rgt = cgmath::vec4((1.0, 0.0, 0.0, 0.0));
    let up = cgmath::vec4((0.0, 1.0, 0.0, 0.0));
    let rotation_axis = cgmath::vec3((0.0, 0.0, -1.0));
    let attitude = CameraAttitude::new(origin, forward, rgt, up, rotation_axis);
    let camera = OrthographicCamera::new(frustum, attitude);

    camera
}

#[test]
fn test_orthographic_camera_z_axis_tans_mat() {
    let camera = orthographic_camera_z_axis();
    let trans_mat = Matrix4::new(
        1.0, 0.0,  0.0, 0.0,
        0.0, 1.0,  0.0, 0.0,
        0.0, 0.0,  1.0, 0.0,
        0.0, 0.0, -5.0, 1.0
    );

    assert_eq!(camera.trans_mat, trans_mat);
}

#[test]
fn test_orthographic_camera_z_axis_rot_mat() {
    let camera = orthographic_camera_z_axis();
    let rot_mat = Matrix4::new(
        1.0, 0.0, 0.0, 0.0,
        0.0, 1.0, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, 1.0
    );

    assert_eq!(camera.rot_mat, rot_mat);
}

#[test]
fn test_orthographic_camera_z_axis_view_mat() {
    let camera = orthographic_camera_z_axis();
    let view_mat = Matrix4::new(
        1.0, 0.0,  0.0,  0.0,
        0.0, 1.0,  0.0,  0.0,
        0.0, 0.0,  1.0,  0.0,
        0.0, 0.0, -5.0,  1.0
    );

    assert_eq!(camera.view_mat, view_mat);
}

#[test]
fn test_orthographic_camera_z_axis_proj_mat() {
    let camera = orthographic_camera_z_axis();
    let proj_mat = Matrix4::new(
        0.25,  0.0,   0.0, 0.0,
        0.0, 0.25,   0.0, 0.0,
        0.0,  0.0, -0.25, 0.0,
        0.0,  0.0,   0.0, 1.0
    );

    assert_eq!(camera.proj_mat, proj_mat);
}
