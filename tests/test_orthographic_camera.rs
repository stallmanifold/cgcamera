use cgcamera::{CameraAttitude, Frustum, OrthographicCamera};
use cgmath::{Vector3, Matrix4};


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


#[test]
fn test_orthographic_camera_should_map_points_inside_frustum_inside_canonical_view_volume() {
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
    let camera = OrthographicCamera::new(frustum, attitude);

    let p_cam = cgmath::vec4((-3.0, 3.0, -3.0, 1.0));
    let p_cvv = camera.proj_mat * p_cam;
    let p_cvv = cgmath::vec3((p_cvv.x, p_cvv.y, p_cvv.z));
    let cvv = AxisAlignedBoundingBox::unit_aabb();

    assert!(cvv.contains(p_cvv), "p_cam = {}; p_cvv = {}", p_cam, p_cvv);
}

#[test]
fn test_orthographic_camera_should_map_points_outside_frustum_outside_canonical_view_volume() {
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
    let camera = OrthographicCamera::new(frustum, attitude);

    let p_cam = cgmath::vec4((300.0, 300.0, 300.0, 1.0));
    let p_cvv = camera.proj_mat * p_cam;
    let p_cvv = cgmath::vec3((p_cvv.x, p_cvv.y, p_cvv.z));
    let cvv = AxisAlignedBoundingBox::unit_aabb();

    assert!(!cvv.contains(p_cvv), "p_cam = {}; p_cvv = {}", p_cam, p_cvv);
}

#[test]
fn test_orthographic_camera_z_axis() {
    let left = -4.0;
    let right = 4.0;
    let bottom = -4.0;
    let top = 4.0;
    let near = 4.0;
    let far = -4.0;
    let frustum = Frustum::new(left, right, bottom, top, near, far);

    let origin = cgmath::vec3((0.0, 0.0, 5.0));
    let look = cgmath::vec4((0.0, 0.0, -1.0, 0.0));
    let rgt = cgmath::vec4((1.0, 0.0, 0.0, 0.0));
    let up = cgmath::vec4((0.0, 1.0, 0.0, 0.0));
    let rotation_axis = cgmath::vec3((0.0, 0.0, -1.0));
    let attitude = CameraAttitude::new(origin, look, rgt, up, rotation_axis);
    let camera = OrthographicCamera::new(frustum, attitude);

    let trans_mat = Matrix4::new(
        1.0, 0.0,  0.0, 0.0,
        0.0, 1.0,  0.0, 0.0,
        0.0, 0.0,  1.0, 0.0,
        0.0, 0.0, -5.0, 1.0
    );

    let view_mat = Matrix4::new(
        1.0, 0.0,  0.0,  0.0,
        0.0, 1.0,  0.0,  0.0,
        0.0, 0.0,  1.0,  0.0,
        0.0, 0.0, -5.0,  1.0
    );

    let proj_mat = Matrix4::new(
        0.25,  0.0,   0.0, 0.0,
         0.0, 0.25,   0.0, 0.0,
         0.0,  0.0, -0.25, 0.0,
         0.0,  0.0,   0.0, 1.0
    );

    assert_eq!(camera.trans_mat, trans_mat);
    assert_eq!(camera.view_mat, view_mat);
    assert_eq!(camera.proj_mat, proj_mat);
}
