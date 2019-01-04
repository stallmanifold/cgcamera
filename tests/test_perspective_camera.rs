use cgmath::{Matrix4, Vector3};
use cgcamera::{
    CameraAttitude, Frustum, FrustumFov, PerspectiveCamera, PerspectiveFovCamera
};


struct AxisAlignedBoundingBox {
    left: f32,
    right: f32,
    bottom: f32,
    top: f32,
    near: f32,
    far: f32,
}

impl AxisAlignedBoundingBox {
    fn new(left: f32, right: f32, bottom: f32, top: f32, near: f32, far: f32) -> Self {
        Self { left, right, bottom, top, near, far }
    }

    fn contains(&self, point: Vector3) -> bool {
        (point.x >= self.left) && (point.x <= self.right) &&
            (point.y >= self.bottom) && (point.y <= self.top) &&
            (point.z >= self.far && point.z <= self.near)
    }

    fn unit_aabb() -> Self {
        Self::new(-1.0, 1.0, -1.0, 1.0, 1.0, -1.0)
    }
}

#[test]
fn test_perspective_camera_model() {
    let left = -4.0;
    let right = 4.0;
    let bottom = -4.0;
    let top = 4.0;
    let near = 0.1;
    let far = 4.0;
    let frustum = Frustum::new(left, right, bottom, top, near, far);

    let origin = cgmath::vec3((0.0, 0.0, 5.0));
    let forward = cgmath::vec4((0.0, 0.0, -1.0, 0.0));
    let right = cgmath::vec4((1.0, 0.0, 0.0, 0.0));
    let up = cgmath::vec4((0.0, 1.0, 0.0, 0.0));
    let rotation_axis = cgmath::vec3((0.0, 0.0, -1.0));
    let attitude = CameraAttitude::new(origin, forward, right, up, rotation_axis);
    let camera = PerspectiveCamera::new(frustum, attitude);

    let p_wor = cgmath::vec4((-3.0, 3.0, 3.0, 1.0));
    let p_cam = camera.view_mat * p_wor;
    let p_cvv = camera.proj_mat * p_cam;
    let p_cvv = cgmath::vec3((p_cvv.x, p_cvv.y, p_cvv.z));
    let cvv = AxisAlignedBoundingBox::unit_aabb();

    assert!(cvv.contains(p_cvv), "p_wor = {}; p_cam = {}; p_cvv = {}", p_wor, p_cam, p_cvv);
}


#[test]
fn test_fov_and_frustum_should_yield_same_camera_model() {
    let origin = cgmath::vec3((0.0, 0.0, 5.0));
    let forward = cgmath::vec4((0.0, 0.0, -1.0, 0.0));
    let right = cgmath::vec4((1.0, 0.0, 0.0, 0.0));
    let up = cgmath::vec4((0.0, 1.0, 0.0, 0.0));
    let rotation_axis = cgmath::vec3((0.0, 0.0, -1.0));
    let attitude = CameraAttitude::new(origin, forward, right, up, rotation_axis);

    let near = 0.1;
    let far = 100.0;
    let fovy = 67.0; // degrees
    let aspect = 720.0 / 480.0;
    let frustum_fov = FrustumFov::new(near, far, fovy, aspect);
    let frustum_box = Frustum::from(frustum_fov);

    let persp_box = PerspectiveCamera::new(frustum_box, attitude);
    let persp_fov = PerspectiveFovCamera::new(frustum_fov, attitude);

    assert_eq!(persp_box.view_mat, persp_fov.view_mat);
    assert_eq!(persp_box.proj_mat, persp_fov.proj_mat);
}

#[test]
fn test_perspective_camera_should_map_points_inside_frustum_inside_canonical_view_volume() {
    let left = -4.0;
    let right = 4.0;
    let bottom = -4.0;
    let top = 4.0;
    let near = 0.1;
    let far = 100.0;
    let frustum = Frustum::new(left, right, bottom, top, near, far);

    let origin = cgmath::vec3((0.0, 0.0, 5.0));
    let forward = cgmath::vec4((0.0, 0.0, -1.0, 0.0));
    let right = cgmath::vec4((1.0, 0.0, 0.0, 0.0));
    let up = cgmath::vec4((0.0, 1.0, 0.0, 0.0));
    let rotation_axis = cgmath::vec3((0.0, 0.0, -1.0));
    let attitude = CameraAttitude::new(origin, forward, right, up, rotation_axis);
    let camera = PerspectiveCamera::new(frustum, attitude);

    let p_cam = cgmath::vec4((-3.0, 3.0, 3.0, 1.0));
    let p_cvv = camera.proj_mat * p_cam;
    let p_cvv = cgmath::vec3((p_cvv.x, p_cvv.y, p_cvv.z));
    let cvv = AxisAlignedBoundingBox::unit_aabb();

    assert!(cvv.contains(p_cvv), "p_cam = {}; p_cvv = {}", p_cam, p_cvv);
}

#[test]
fn test_perspective_camera_should_map_points_outside_frustum_outside_canonical_view_volume() {
    let left = -4.0;
    let right = 4.0;
    let bottom = -4.0;
    let top = 4.0;
    let near = 0.1;
    let far = 4.0;
    let frustum = Frustum::new(left, right, bottom, top, near, far);

    let origin = cgmath::vec3((0.0, 0.0, 5.0));
    let forward = cgmath::vec4((0.0, 0.0, -1.0, 0.0));
    let right = cgmath::vec4((1.0, 0.0, 0.0, 0.0));
    let up = cgmath::vec4((0.0, 1.0, 0.0, 0.0));
    let rotation_axis = cgmath::vec3((0.0, 0.0, -1.0));
    let attitude = CameraAttitude::new(origin, forward, right, up, rotation_axis);
    let camera = PerspectiveCamera::new(frustum, attitude);

    let p_cam = cgmath::vec4((300.0, 300.0, 300.0, 1.0));
    let p_cvv = camera.proj_mat * p_cam;
    let p_cvv = cgmath::vec3((p_cvv.x, p_cvv.y, p_cvv.z));
    let cvv = AxisAlignedBoundingBox::unit_aabb();

    assert!(!cvv.contains(p_cvv), "p_cam = {}; p_cvv = {}", p_cam, p_cvv);
}

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
        1.0,   0.0,  0.0,  0.0,
        0.0,   1.0,  0.0,  0.0,
        0.0,   0.0,  0.0, -4.0,
        0.0,   0.0,  1.0,  0.0
    );

    assert_eq!(camera.proj_mat, proj_mat);
}
