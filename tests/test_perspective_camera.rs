use cgmath::{Matrix4, Vector3};
use cgcamera::{
    CameraAttitude, Frustum, FrustumFov, PerspectiveCamera, PerspectiveFovCamera
};


fn perspective_camera_z_axis() -> PerspectiveCamera {
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
    let camera = PerspectiveCamera::new(frustum, attitude);

    camera
}

#[test]
fn test_perspective_camera_z_axis_trans_mat() {
    let camera = perspective_camera_z_axis();
    let trans_mat = Matrix4::new(
        1.0, 0.0,   0.0, 0.0,
        0.0, 1.0,   0.0, 0.0,
        0.0, 0.0,   1.0, 0.0,
        0.0, 0.0,  -5.0, 1.0
    );

    assert_eq!(camera.trans_mat, trans_mat);
}

#[test]
fn test_perspective_camera_z_axis_rot_mat() {
    let camera = perspective_camera_z_axis();
    let rot_mat = Matrix4::new(
        -1.0,  0.0, 0.0, 0.0,
         0.0, -1.0, 0.0, 0.0,
         0.0,  0.0, 1.0, 0.0,
         0.0,  0.0, 0.0, 1.0
    );

    assert_eq!(camera.rot_mat, rot_mat);
}

#[test]
fn test_perspective_camera_z_axis_view_mat() {
    let camera = perspective_camera_z_axis();
    let view_mat = Matrix4::new(
        -1.0,  0.0,  0.0,  0.0,
         0.0, -1.0,  0.0,  0.0,
         0.0,  0.0,  1.0,  0.0,
         0.0,  0.0, -5.0,  1.0
    );

    assert_eq!(camera.view_mat, view_mat);
}

#[test]
fn test_perspective_camera_z_axis_proj_mat() {
    let camera = perspective_camera_z_axis();
    let proj_mat = Matrix4::new(
        1.0,  0.0,  0.0,  0.0,
        0.0,  1.0,  0.0,  0.0,
        0.0,  0.0,  0.0,  1.0,
        0.0,  0.0, -4.0,  0.0
    );

    assert_eq!(camera.proj_mat, proj_mat);
}

fn perspective_camera_x_axis() -> PerspectiveCamera {
    let left = -4.0;
    let right = 4.0;
    let bottom = -4.0;
    let top = 4.0;
    let near = 4.0;
    let far = -4.0;
    let frustum = Frustum::new(left, right, bottom, top, near, far);

    let origin = cgmath::vec3((5.0, 0.0, 0.0));
    let forward = cgmath::vec4((1.0, 0.0, 0.0, 0.0));
    let rgt = cgmath::vec4((0.0, 1.0, 0.0, 0.0));
    let up = cgmath::vec4((0.0, 0.0, 1.0, 0.0));
    let rotation_axis = cgmath::vec3((-1.0, 0.0, 0.0));
    let attitude = CameraAttitude::new(origin, forward, rgt, up, rotation_axis);
    let camera = PerspectiveCamera::new(frustum, attitude);

    camera
}

#[test]
fn test_perspective_camera_x_axis_trans_mat() {
    let camera = perspective_camera_x_axis();
    let trans_mat = Matrix4::new(
         1.0, 0.0,  0.0, 0.0,
         0.0, 1.0,  0.0, 0.0,
         0.0, 0.0,  1.0, 0.0,
        -5.0, 0.0,  0.0, 1.0
    );

    assert_eq!(camera.trans_mat, trans_mat);
}

#[test]
fn test_perspective_camera_x_axis_rot_mat() {
    let camera = perspective_camera_x_axis();
    let rot_mat = Matrix4::new(
        1.0,  0.0,  0.0, 0.0,
        0.0, -1.0,  0.0, 0.0,
        0.0,  0.0, -1.0, 0.0,
        0.0,  0.0,  0.0, 1.0
    );

    assert_eq!(camera.rot_mat, rot_mat);
}

#[test]
fn test_perspective_camera_x_axis_view_mat() {
    let camera = perspective_camera_x_axis();
    let view_mat = Matrix4::new(
         1.0,  0.0,  0.0,  0.0,
         0.0, -1.0,  0.0,  0.0,
         0.0,  0.0, -1.0,  0.0,
        -5.0,  0.0,  0.0,  1.0
    );

    assert_eq!(camera.view_mat, view_mat);
}

#[test]
fn test_perspective_camera_x_axis_proj_mat() {
    let camera = perspective_camera_x_axis();
    let proj_mat = Matrix4::new(
        1.0,  0.0,  0.0,  0.0,
        0.0,  1.0,  0.0,  0.0,
        0.0,  0.0,  0.0,  1.0,
        0.0,  0.0, -4.0,  0.0
    );

    assert_eq!(camera.proj_mat, proj_mat);
}
