/**
 * Copyright (C) 2015 by Liangliang Nan (liangliang.nan@gmail.com)
 * https://3d.bk.tudelft.nl/liangliang/
 *
 * This file is part of Easy3D. If it is useful in your research/work,
 * I would be grateful if you show your appreciation by citing it:
 * ------------------------------------------------------------------
 *      Liangliang Nan.
 *      Easy3D: a lightweight, easy-to-use, and efficient C++
 *      library for processing and rendering 3D data. 2018.
 * ------------------------------------------------------------------
 * Easy3D is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License Version 3
 * as published by the Free Software Foundation.
 *
 * Easy3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "triangulation.h"
#include "matrix_algo.h"
#include <easy3d/optimizer/optimizer_lm.h>
#include <chrono>
// added tuple class to be able to return normalized points and their T matrices at the same time
#include <tuple>


using namespace easy3d;

// given:

/// convert a 3 by 3 matrix of type 'Matrix<double>' to mat3
mat3 to_mat3(Matrix<double> &M) {
    mat3 result;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j)
            result(i, j) = M(i, j);
    }
    return result;
}


/// convert M of type 'matN' (N can be any positive integer) to type 'Matrix<double>'
template<typename mat>
Matrix<double> to_Matrix(const mat &M) {
    const int num_rows = M.num_rows();
    const int num_cols = M.num_columns();
    Matrix<double> result(num_rows, num_cols);
    for (int i = 0; i < num_rows; ++i) {
        for (int j = 0; j < num_cols; ++j)
            result(i, j) = M(i, j);
    }
    return result;
}

// own functions:

/// check input validity
void test_input(const std::vector<vec3> &points_1, const std::vector<vec3> &points_2){



    std::vector<int> duplicateLocations;

    for (unsigned int i = 0; i < points_1.size(); i ++ ) {
        // check for negative 3D points
        if (points_1[i].x < 0 || points_1[i].y < 0  ||
            points_2[i].x < 0 || points_2[i].y < 0 ) {
            std::cout << "Invalid 2d point with negative coordinates: ("
                      << points_1[i].x << " "
                      << points_1[i].y << ") at location: "
                      << i << std::endl;

            duplicateLocations.emplace_back(i);
            continue;
        }

        // check if points are homogeneous with z = 1.0
        if (points_1[i].z != 1 || points_2[i].z != 1){
            std::cout << "Invalid constructed 2d point: ("
                      << points_1[i].x << " "
                      << points_1[i].z << ") at location: "
                      << i << std::endl;

            duplicateLocations.emplace_back(i);
            continue;
        }

        // check for duplicates
        for (unsigned int j = 0; j < points_1.size(); j ++) {
            if ( i >= j) {
                continue;
            }
            if (points_1[i].x == points_1[j].x
                && points_1[i].y == points_1[j].y) {

                std::cout << "Duplicate 2d point point with coordinates: ("
                          << points_1[i][0] << " "
                          << points_1[i][1] << ") - ("
                          << points_2[i][0] << " "
                          << points_2[i][1] << ") at location: "
                          << i << std::endl;

                duplicateLocations.emplace_back(i);
                break;
            }
        }
    }


    // check if the size after removal is big enough
    if (points_1.size() - duplicateLocations.size() < 8){
        std::cout << "=========================================================" << std::endl;
        std::cout << "expecting at least 8 unique pairs of corresponding points" << std::endl;
        std::cout << "=========================================================" << std::endl;
    }
}


/// Normalize input points
std::tuple<std::vector<vec3>, std::vector<vec3>, mat3, mat3> normalize_input_all(const std::vector<vec3> &points0,
                                                                                 const std::vector<vec3> &points1){
    std::vector<vec3> pts_norm0;
    std::vector<vec3> pts_norm1;

    /// step 1: translation
    // the origin of the new coordinate system should be located at the centroid of the image points

    double x_coords0 = 0;    // total x coordinates
    double y_coords0 = 0;    // total y coordinates
    double z_coords0 = 0;    // total z coordinates

    double x_coords1 = 0;    // total x coordinates
    double y_coords1 = 0;    // total y coordinates
    double z_coords1 = 0;    // total z coordinates

    for (int i = 0; i < points0.size(); ++i) {
        // average x, y, (z) = centroid
        x_coords0 += points0[i].x;
        x_coords1 += points1[i].x;

        y_coords0 += points0[i].y;
        y_coords1 += points1[i].y;

        z_coords0 += points0[i].z;
        z_coords1 += points1[i].z;
    }

    double x_av0 = x_coords0 / (points0.size());     // centroid x
    double x_av1 = x_coords1 / (points0.size());     // centroid x
    double y_av0 = y_coords0 / (points0.size());     // centroid y
    double y_av1 = y_coords1 / (points0.size());     // centroid y
    double z_av0 = z_coords0 / (points0.size());     // centroid z
    double z_av1 = z_coords1 / (points0.size());     // centroid z

    // total distance to origin (centroid), after translation, per image
    // image 0
    double dist0 = 0;
    for (vec3 point : points0){
        double x_diff = point.x - x_av0;
        double y_diff = point.y - y_av0;
        double z_diff = point.z - z_av0;

        dist0 += sqrt(pow(x_diff, 2) + pow(y_diff, 2) + pow(z_diff, 2));
    }
    // image 1
    double dist1 = 0;
    for (vec3 point : points1){
        double x_diff = point.x - x_av1;
        double y_diff = point.y - y_av1;
        double z_diff = point.z - z_av1;

        dist1 += sqrt(pow(x_diff, 2) + pow(y_diff, 2) + pow(z_diff, 2));
    }

    // average distance to centroid, per image
    double dist0_av = dist0 / points0.size();
    double dist1_av = dist1 / points0.size();

    /// step 2: scaling and translating
    mat3 T0(1.0f);           // translation transformation matrix, diagonal = 1 & rest = 0
    mat3 T1(1.0f);

    mat3 T0_scale(1.0f);
    mat3 T1_scale(1.0f);

    mat3 T0_trans(1.0f);
    mat3 T1_trans(1.0f);

    // the mean square distance of the transformed image points from the origin should be 2 pixels
    double scaling_factor0 = sqrt(2) / dist0_av;
    double scaling_factor1 = sqrt(2) / dist1_av;

    T0_trans(0, 2) = - x_av0;   // x translation
    T0_trans(1, 2) = - y_av0;   // y translation

    T1_trans(0, 2) = - x_av1;   // x translation
    T1_trans(1, 2) = - y_av1;   // y translation

    T0_scale(0,0) = scaling_factor0;  // x scaling
    T0_scale(1,1) = scaling_factor0;  // y scaling

    T1_scale(0,0) = scaling_factor1;  // x scaling
    T1_scale(1,1) = scaling_factor1;  // y scaling
    // z scaling/translation should just remain 1: homogeneous coordinates

    T0 = T0_scale * T0_trans;
    T1 = T1_scale * T1_trans;

    /// normalize the points
    // apply transformations to all points
    for (int i = 0; i < points0.size(); ++i) {
        vec3 res0 = T0 * points0[i];
        pts_norm0.push_back(res0 / res0.z);

        vec3 res1 = T1 * points1[i];
        pts_norm1.push_back(res1 / res1.z);
    }

    std::tuple<std::vector<vec3>, std::vector<vec3>, mat3, mat3> res = make_tuple(pts_norm0, pts_norm1, T0, T1);
    return res;
};


/// Step 1: estimate fundamental matrix F
mat3 estimate_F(std::vector<vec3> pts0_norm, std::vector<vec3> pts1_norm){
    /// construct linear system W
    // W size = n * 9
    Matrix<double> W(pts0_norm.size(), 9, 0.0);

    // fill W
    for (int i = 0; i < pts0_norm.size(); ++i) {
        vec3 pt0 = pts0_norm[i];
        vec3 pt1 = pts1_norm[i];
        std::vector<double> elem(9, 1.0);

        elem[0] = pt0.x * pt1.x;
        elem[1] = pt0.y * pt1.x;
        elem[2] = pt1.x;
        elem[3] = pt0.x * pt1.y;
        elem[4] = pt0.y * pt1.y;
        elem[5] = pt1.y;
        elem[6] = pt0.x;
        elem[7] = pt0.y;

        W.set_row(elem, i);
    }

    /// solve Wf = 0 --> F'
    int m = pts0_norm.size();
    int n = 9;

    Matrix<double> U_W(m, m, 0.0);   // initialized with 0s
    Matrix<double> S_W(m, n, 0.0);
    Matrix<double> V_W(n, n, 0.0);

    // SVD decomposition of W
    svd_decompose(W, U_W, S_W, V_W);

    // F' from f (f = last column of V)
    Matrix<double> F_e(3, 3, 0.0);

    Matrix<double> V_W_T = transpose(V_W);
    // populate F'
    for (int i = 0; i < 3; i ++){
        for (int j = 0; j < 3; j++){
            F_e(i,j) = V_W(i * 3 + j, 8);
        }
    }

    /// minimization & rank 2 --> F

    // decompose F' (3 x 3)
    Matrix<double> U_F(3, 3, 0.0);   // initialized with 0s
    Matrix<double> S_F(3, 3, 0.0);   // initialized with 0s
    Matrix<double> V_F(3, 3, 0.0);   // initialized with 0s

    // SVD decomposition of F'
    svd_decompose(F_e, U_F, S_F, V_F);

    // enforce rank 2
    S_F(2, 2) = 0.0;

    /// construct final F
    Matrix<double> svd_res = U_F * S_F * transpose(V_F);
    mat3 F = to_mat3(svd_res);  // to fixed size for efficiency

    return F;
}


/// Step 3: Determine the 3D coordinates
// linear triangulation method
vec3 triangulate_linear(vec3 point0, vec3 point1, mat3 K, mat3 R, vec3 t){

    /// projection matrices of the two cameras
    // M & M', named M0 & M1 here

    // M = K[I 0]
    mat34 M0(1.0f);
    M0(2,3) = 0.0;
    M0 = K * M0;

    // M' = K'[R t]
    mat34 M1(1.0f);
    for (int i = 0; i < R.num_rows(); ++i) {
        vec3 row = R.row(i);
        M1.set_row(i, vec4(row[0], row[1], row[2], 0.0));
    }
    M1.set_col(3, t);
    M1 = K * M1;

    /// construct matrix A
    // AP = 0
    Matrix<double> A(4, 4, 0.0);
    vec4 elem00 = point0.x * M0.row(2) - M0.row(0);
    vec4 elem10 = point0.y * M0.row(2) - M0.row(1);
    vec4 elem20 = point1.x * M1.row(2) - M1.row(0);
    vec4 elem30 = point1.y * M1.row(2) - M1.row(1);

    // fill A
    for (int i = 0; i < 4; ++i) {
        A(0,i) = elem00[i];
        A(1,i) = elem10[i];
        A(2,i) = elem20[i];
        A(3,i) = elem30[i];
    }

    // SVD of A
    int h = A.rows();
    int w = A.cols();

    Matrix<double> U(h, h, 0.0);
    Matrix<double> S(h, w, 0.0);
    Matrix<double> V(w, w, 0.0);

    svd_decompose(A, U, S, V);

    /// construct point P
    Matrix<double> P(4, 1, 0.0);
    for (int i = 0; i < 4; ++i) {
        P(i,0) = V(i, w-1);
    }

    // normalize P (homogeneous coordinates, w = 1.0)
    Matrix<double> res = (P / P(3,0));

    // 4D (hom) to 3D (cartesian)
    vec3 point_3d = {float(res(0,0)), float(res(1,0)), float(res(2,0))};

    return point_3d;
}

/// non-linear triangulation with Gauss-Newton method
// residual error calculation
vec3 residual(vec3 point, mat34 M, vec4 point_3d){

    //reprojected point
    vec3 p0_rep_hom = M * point_3d;

    // homogenous with z=1.0
    vec3 prep = p0_rep_hom / p0_rep_hom.z;

    vec3 error = prep - point;

    return error;
};

// sum of squared errors (min (|| r(x) ||^2))
double sum_square_error(vec3 point0, vec3 point1, mat34 M0, mat34 M1, vec4 P_est_curr){
    vec3 error0 = residual(point0, M0, P_est_curr);
    vec3 error1 = residual(point1, M1, P_est_curr);

    double total_error = pow((norm(error0) + norm(error1)), 2);

    return total_error;
};


// non-linear triangulation method
vec3 triangulate_nonlinear(vec3 point0, vec3 point1, vec3 point_3d, mat3 K, mat3 R, vec3 t){
    // p = M P               p' = M' P

    // minimize P_est : dist(M P_est - p)^2 + dist(M' P_est - p')^2
    /// projection matrices of the two cameras
    // M & M', named M0 & M1 here

    // M = K[I 0]
    mat34 M0(1.0f);
    M0(2,3) = 0.0;
    M0 = K * M0;

    // M' = K'[R t]
    mat34 M1(1.0f);
    for (int i = 0; i < R.num_rows(); ++i) {
        vec3 row = R.row(i);
        M1.set_row(i, vec4(row[0], row[1], row[2], 0.0));
    }
    M1.set_col(3, t);
    M1 = K * M1;

    // 3d point in homogeneous coordinates
    vec4 P_est = {point_3d.x, point_3d.y, point_3d.z, 1.0};

    /// construct residuals e
    // residuals per point
    vec3 e0 = residual(point0, M0, P_est);
    vec3 e1 = residual(point1, M1, P_est);

    // residual matrix e (4x1)
    Matrix<double> e(4,1);
    e.set_column({e0.x, e0.y, e1.x, e1.y}, 0);

    /// construct Jacobian J (4x3)
    Matrix<double> J(4,3);

    // delta P coor is always 1.0
    vec4 Px_delta = {P_est.x+1, P_est.y, P_est.z, P_est.w};
    vec4 Py_delta = {P_est.x, P_est.y+1, P_est.z, P_est.w};
    vec4 Pz_delta = {P_est.x, P_est.y, P_est.z+1, P_est.w};

    vec3 e0_delta_px = residual(point0, M0, Px_delta);
    vec3 e0_delta_py = residual(point0, M0, Py_delta);
    vec3 e0_delta_pz = residual(point0, M0, Pz_delta);

    vec3 e1_delta_px = residual(point1, M1, Px_delta);
    vec3 e1_delta_py = residual(point1, M1, Py_delta);
    vec3 e1_delta_pz = residual(point1, M1, Pz_delta);

    // point 0

    for (int i = 0; i < 2; ++i) {
        J (i, 0) = e0_delta_px[i];
        J (i, 1) = e0_delta_py[i];
        J (i, 2) = e0_delta_pz[i];
    }

    for (int i = 0; i < 2; ++i) {
        J (i + 2, 0) = e1_delta_px[i];
        J (i + 2, 1) = e1_delta_py[i];
        J (i + 2, 2) = e1_delta_pz[i];
    }

    /// determine delta(P)
    Matrix<double> JTJ = transpose(J) * J;
    Matrix<double> JTJ_inv;
    inverse(JTJ, JTJ_inv);
    Matrix<double> delta_P = - JTJ_inv * transpose(J) * e;
    vec4 delta_P_vec = {float(delta_P(0,0)), float(delta_P(1,0)), float(delta_P(2,0)), 0};

    /// iterate until error is acceptable, or iteration was done x times
    int steps = 50;             // total number of iterations
    double epsilon = 0.00001;   // maximal acceptable error

    int step_curr = 0;

    double error_curr = sum_square_error(point0, point1, M0, M1, P_est);
    double error_min = error_curr;

    bool cont = true;

    while(cont){
        step_curr ++;
        // take a step with delta p
        vec4 P_est_curr = P_est + delta_P_vec;
        error_curr = sum_square_error(point0, point1, M0, M1, P_est_curr);

        if (error_curr < error_min){
            P_est = P_est_curr;
            error_min = error_curr;
        }
        if (step_curr >= steps){
            cont = false;
        }
        if (error_min <= epsilon){
            cont = false;
        }
    }

    /// homogeneous to cartesian coordinates
    vec4 res = P_est / P_est.w;
    vec3 point_3d_est = {res.x, res.y, res.z};
    return point_3d_est;
}


/// Step 2: Recover relative pose (R & t)
std::tuple<mat3, vec3> correct_direction_poses(mat3 R01, mat3 R02, vec3 u3, mat3 K,
                                               std::vector<vec3> pts0,
                                               std::vector<vec3> pts1){

    // determinant of R has to be positive
    mat3 R1 = determinant(R01) * R01;
    mat3 R2 = determinant(R02) * R02;

    // possible ts
    vec3 t1 = u3;
    vec3 t2 = -u3;

    // for all 4 possible combinations of R and t
    // see which one gives the most 3d points that are in front of both cameras
    mat3 R_best;
    vec3 t_best;

    int cnt_max = 0;
    mat3 R_curr = R1;
    vec3 t_curr = t1;

    for (int sit = 0; sit < 4; ++sit) {
        int cnt = 0;

        // situation selection:
        // sit 0 = R1 t1
        // sit 1 = R2 t1
        // sit 2 = R1 t2
        // sit 3 = R2 t2
        if (sit == 1){
            R_curr = R2;
        }
        if (sit == 2){
            R_curr = R1;
            t_curr = t2;
        }
        if (sit == 3){
            R_curr = R2;
        }
        for (int i = 0; i < pts0.size(); ++i) {
            // triangulate the point-pair
            vec3 pt_3d = triangulate_linear(pts0[i], pts1[i], K, R_curr, t_curr);

            // find if the point is in front of both cameras
            vec3 pt_3d_1 = R_curr * pt_3d + t_curr;     // camera center of 2nd camera
            if (pt_3d.z > 0 && pt_3d_1.z > 0){
                cnt ++;
            }
        }
        // check if more points are in front of both cameras here
        if (cnt > cnt_max){
            cnt_max = cnt;
            R_best = R_curr;
            t_best = t_curr;
        }
    }
    return std::make_tuple(R_best, t_best);
};


std::tuple<mat3, vec3> relative_position(mat3 E, mat3 K, std::vector<vec3> points_0, std::vector<vec3> points_1){
    /// determine relative position image 1 to image 0
    // as translation t and rotation R

    // W =  |   0  -1   0   |      Z =  |   0   1   0   |
    //      |   1   0   0   |           |  -1   0   0   |
    //      |   0   0   1   |           |   0   0   0   |

    mat3 W = {0, -1, 0, 1, 0, 0, 0, 0, 1};

    // E = U S V^T          [t]x = U Z U^T
    // E = [t]x R           t = +- u3

    // R = U W V^T or U W^T V^T

    /// SVD decomposition of E
    Matrix<double> U_mat(3, 3, 0.0);
    Matrix<double> S_mat(3, 3, 0.0);
    Matrix<double> V_mat(3, 3, 0.0);

    Matrix<double> E_mat = to_Matrix(E);
    svd_decompose(E_mat, U_mat, S_mat, V_mat);

    mat3 U = to_mat3(U_mat);
    mat3 V = to_mat3(V_mat);

    /// rotation R
    mat3 R1 = U * W * transpose(V);
    mat3 R2 = U * transpose(W) * transpose(V);

    /// translation t
    vec3 u3 = U.col(2);

    /// correct direction
    // of both R and t
    std::tuple<mat3, vec3> rel_pos_corr = correct_direction_poses(R1, R2, u3, K, points_0, points_1);

    return rel_pos_corr;
};



bool Triangulation::triangulation(
        float fx, float fy,     /// input: the focal lengths (same for both cameras)
        float cx, float cy,     /// input: the principal point (same for both cameras)
        const std::vector<vec3> &points_0,    /// input: image points (in homogenous coordinates) in the 1st image.
        const std::vector<vec3> &points_1,    /// input: image points (in homogenous coordinates) in the 2nd image.
        std::vector<vec3> &points_3d,         /// output: reconstructed 3D points
        mat3 &R,   /// output: recovered rotation of 2nd camera (used for updating the viewer and visual inspection)
        vec3 &t    /// output: recovered translation of 2nd camera (used for updating the viewer and visual inspection)
) const
{
    //--------------------------------------------------------------------------------------------------------------
    // implementation starts ...

    auto startTime = std::chrono::high_resolution_clock::now();

    /// check if the input is valid (always good because you never known how others will call your function).
    test_input(points_0, points_1);

    /// estimate the fundamental matrix F
    // normalize input points
    std::tuple<std::vector<vec3>, std::vector<vec3>, mat3, mat3> norm = normalize_input_all(points_0, points_1);
    std::vector<vec3> pts0_norm = std::get<0>(norm);
    std::vector<vec3> pts1_norm = std::get<1>(norm);
    mat3 T_0 = std::get<2>(norm);
    mat3 T_1 = std::get<3>(norm);

    // construct F
    mat3 F_init = estimate_F(pts0_norm, pts1_norm);
    // de-normalize F
    mat3 F_norm = transpose(T_1) * F_init * T_0;
    // rescale F
    mat3 F = F_norm / F_norm(2,2);

    /// construct intrinsic matrix K from given input parameters
    // K =  |   fx    skew    cx   |
    //      |   0      fy     cy   |
    //      |   0      0       1   |
    mat3 K(1.0f);
    K(0,0) = fx;
    K(1,1) = fy;
    K(0,2) = cx;
    K(1,2) = cy;

    /// compute the essential matrix E;
    // E = K^T F K  --> intrinsics are known
    // K' = K according to assignment
    mat3 E = transpose(K) * F * K;

    /// recover rotation R and t.
    // relative position between the two cameras
    std::tuple<mat3, vec3> pos = relative_position(E, K, points_0, points_1);
    R = std::get<0>(pos);
    t = std::get<1>(pos);

    /// Reconstruct 3D points. The main task is triangulation
    // triangulate a pair of image points
    // (i.e., compute the 3D coordinates for each corresponding point pair)
    for (int i = 0; i < points_0.size(); ++i) {
        vec3 pt3d = triangulate_linear(points_0[i], points_1[i], K, R, t);

        // non-linear improvement on the linear triangulation
        vec3 pt3d_imp = triangulate_nonlinear(points_0[i], points_1[i], pt3d, K, R, t);

        points_3d.push_back(pt3d_imp);
    }

    auto endTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsedTime = endTime - startTime;
    std::cout << "triagulation duration: " << elapsedTime.count() << " s" << std::endl;

    // return true or false on success or failure
    return points_3d.size() > 0;
}