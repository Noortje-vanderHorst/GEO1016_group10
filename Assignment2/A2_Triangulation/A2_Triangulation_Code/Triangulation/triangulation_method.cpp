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

// todo: added this, don't know if it's allowed though
#include <tuple>


using namespace easy3d;


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


/// comments and hints provided (we can just throw this whole thing away in the final version)
void Comments(){
    /// NOTE: there might be multiple workflows for reconstructing 3D geometry from corresponding image points.
    ///       This assignment uses the commonly used one explained in our lecture.
    ///       It is advised to define a function for each sub-task. This way you have a clean and well-structured
    ///       implementation, which also makes testing and debugging easier. You can put your other functions above
    ///       triangulation(), or feel free to put them in one or multiple separate files.

    std::cout << "\nTODO: I am going to implement the triangulation() function in the following file:" << std::endl
              << "\t    - triangulation_method.cpp\n\n";

    std::cout << "[Liangliang]:\n"
                 "\tFeel free to use any data structure and function offered by Easy3D, in particular the following two\n"
                 "\tfiles for vectors and matrices:\n"
                 "\t    - easy3d/core/mat.h  Fixed-size matrices and related functions.\n"
                 "\t    - easy3d/core/vec.h  Fixed-size vectors and related functions.\n"
                 "\tFor matrices with unknown sizes (e.g., when handling an unknown number of corresponding points\n"
                 "\tstored in a file, where their sizes can only be known at run time), a dynamic-sized matrix data\n"
                 "\tstructure is necessary. In this case, you can use the templated 'Matrix' class defined in\n"
                 "\t    - Triangulation/matrix.h  Matrices of arbitrary dimensions and related functions.\n"
                 "\tPlease refer to the corresponding header files for more details of these data structures.\n\n"
                 "\tIf you choose to implement the non-linear method for triangulation (optional task). Please refer to\n"
                 "\t'Tutorial_NonlinearLeastSquares/main.cpp' for an example and some explanations. \n\n"
                 "\tIn your final submission, please\n"
                 "\t    - delete ALL unrelated test or debug code and avoid unnecessary output.\n"
                 "\t    - include all the source code (original code framework + your implementation).\n"
                 "\t    - do NOT include the 'build' directory (which contains the intermediate files in a build step).\n"
                 "\t    - make sure your code compiles and can reproduce your results without any modification.\n\n" << std::flush;

    /// Easy3D provides fixed-size matrix types, e.g., mat2 (2x2), mat3 (3x3), mat4 (4x4), mat34 (3x4).
    /// To use these matrices, their sizes should be known to you at the compile-time (i.e., when compiling your code).
    /// Once defined, their sizes can NOT be changed.
    /// In 'Triangulation/matrix.h', another templated 'Matrix' type is also provided. This type can have arbitrary
    /// dimensions and their sizes can be specified at run-time (i.e., when executing your program).
    /// Below are a few examples showing some of these data structures and related APIs.

    /// ----------- fixed-size matrices

    /// define a 3 by 4 matrix M (you can also define 3 by 4 matrix similarly)
    mat34 M(1.0f);  /// entries on the diagonal are initialized to be 1 and others to be 0.

    /// set the first row of M
    M.set_row(0, vec4(1,1,1,1));    /// vec4 is a 4D vector.

    /// set the second column of M
    M.set_col(1, vec4(2,2,2,2));

    /// get the 3 rows of M
    vec4 M1 = M.row(0);
    vec4 M2 = M.row(1);
    vec4 M3 = M.row(2);

    /// ----------- fixed-size vectors

    /// how to quickly initialize a std::vector
    std::vector<double> rows = {0, 1, 2, 3,
                                4, 5, 6, 7,
                                8, 9, 10, 11};
    /// get the '2'-th row of M
    const vec4 b = M.row(2);    // it assigns the requested row to a new vector b

    /// get the '1'-th column of M
    const vec3 c = M.col(1);    // it assigns the requested column to a new vector c

    /// modify the element value at row 2 and column 1 (Note the 0-based indices)
    M(2, 1) = b.x;

    /// apply transformation M on a 3D point p (p is a 3D vector)
    vec3 p(222, 444, 333);
    vec3 proj = M * vec4(p, 1.0f);  // use the homogenous coordinates. result is a 3D vector

    /// the length of a vector
    float len = p.length();
    /// the squared length of a vector
    float sqr_len = p.length2();

    /// the dot product of two vectors
    float dot_prod = dot(p, proj);

    /// the cross product of two vectors
    vec3 cross_prod = cross(p, proj);

    /// normalize this vector
    cross_prod.normalize();

    /// a 3 by 3 matrix (all entries are intentionally NOT initialized for efficiency reasons)
    mat3 F;
    /// ... here you compute or initialize F.
    /// compute the inverse of K
    mat3 invF = inverse(F);

    /// ----------- dynamic-size matrices

    /// define a non-fixed size matrix
    Matrix<double> W(2, 3, 0.0); // all entries initialized to 0.0.

    /// set its first row by a 3D vector (1.1, 2.2, 3.3)
    W.set_row({ 1.1, 2.2, 3.3 }, 0);   // here "{ 1.1, 2.2, 3.3 }" is of type 'std::vector<double>'

    /// get the last column of a matrix
    std::vector<double> last_column = W.get_column(W.cols() - 1);

    // TODO: delete all above demo code in the final submission
}


/// check input validity
void test_input(){};


/// Normalize input points
std::tuple<std::vector<vec3>,  mat3> normalize_input(const std::vector<vec3> &points){
    std::vector<vec3> pts_norm;  // same size as input

    /// step 1: translation
    // the origin of the new coordinate system should be located at the centroid of the image points

    double x_coords = 0;    // total x coordinates
    double y_coords = 0;    // total y coordinates
    // not for z_coords, in this case the average should just be 1

    for (vec3 point : points){
        // average x, y, (z) = centroid
        x_coords += point.x;
        y_coords += point.y;
    }

    double x_av = x_coords / points.size();     // centroid x
    double y_av = y_coords / points.size();     // centroid y

    double dist = 0;        // total distance to origin (centroid), after translation
    for (vec3 point : points){
        double x_diff = point.x - x_av;
        double y_diff = point.y - y_av;

        dist += sqrt(pow(x_diff, 2) + pow(y_diff, 2));
    }

    double dist_av = dist / points.size();

    /// step 2: scaling and translating
    mat3 T(1.0f);           // translation transformation matrix, diagonal = 1 & rest = 0

    mat3 T_scale(1.0f);
    mat3 T_trans(1.0f);

    // the mean square distance of the transformed image points from the origin should be 2 pixels
    double scaling_factor = sqrt(2) / dist_av;

    T_trans(0, 2) = - x_av;   // x translation
    T_trans(1, 2) = - y_av;   // y translation
//    T_trans(2, 2) =  -1.0;   // z translation

    T_scale(0,0) = scaling_factor;  // x scaling
    T_scale(1,1) = scaling_factor;  // y scaling
    T_scale(2,2) = scaling_factor;  // z scaling

    // z scaling/translation should just remain 1
    T = T_scale * T_trans;

    /// normalize the points
    // apply transformations to all points
    for (vec3 point : points){
//        std::cout << "point norm:\n" << point << std::endl;
//        std::cout << "point norm * T:\n" << T * point << std::endl;

//        Matrix<double> T_mat = to_Matrix(T);
//        Matrix<double> point_mat(3,1);
//        point_mat.set_column({point.x, point.y, point.z}, 0);
//        Matrix<double> res_hom = T_mat * point_mat;
//        vec3 res = {float(res_hom(0,0)), float(res_hom(1,0)), float(res_hom(2,0))};

        vec3 res = T * point;
//        std::cout << "not hom:\n" << res_hom / res_hom.z << std::endl;
//        std::cout << "point res:\n" << res << std::endl;

        // todo: hom to cartessian according to liangliang
        pts_norm.push_back(res /  res.z);
    }

    std::tuple<std::vector<vec3>, mat3> res = make_tuple(pts_norm, T);
    return res;
};


/// Step 1: estimate fundamental matrix F
mat3 estimate_F(std::vector<vec3> pts0_norm, std::vector<vec3> pts1_norm){

    /// construct linear system W
    // row of W:
    // [ u u'       v u'        u'      u v'       v v'       v'      u       v       1   ]
    // [ x0 x1      y0 x1       x1      x0 y1      y0 y1      y1      x0      y0      1   ]

    //      0         1         2         3         4         5       6       7       8

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
    // decompose matrix of size (m x n):
    // U = (m x m)
    // S = (m x n)
    // V = (n x n)

    int m = pts0_norm.size();
    int n = 9;

    Matrix<double> U_W(m, m, 0.0);   // initialized with 0s
    Matrix<double> S_W(m, n, 0.0);   // initialized with 0s
    Matrix<double> V_W(n, n, 0.0);   // initialized with 0s

    // SVD decomposition of W
    svd_decompose(W, U_W, S_W, V_W);

    // F' from f (f = last column of V)
    Matrix<double> F_e(3, 3, 0.0);    // initialized with 0s

    // populate F'
    for (int i = 0; i < 3; i ++){
        for (int j = 0; j < 3; j++){
            F_e(i,j) = V_W(i * 3 + j, 8);
        }
    }

//    std::cout << "V_W:\n" << V_W << std::endl;
//    std::cout << "F': " << F_e << std::endl;

    /// minimization & rank 2 --> F
    // F' = U S' V^T
    // where F' is minimized to approach F, and S has rank 2:
    // S' --> S:
    // S =  | s1   0    0  |
    //      | 0    s2   0  |
    //      | 0    0   [0] | <-- s3 = 0
    // F = U S V^T

    // decompose F' (3 x 3)
    Matrix<double> U_F(3, 3, 0.0);   // initialized with 0s
    Matrix<double> S_F(3, 3, 0.0);   // initialized with 0s
    Matrix<double> V_F(3, 3, 0.0);   // initialized with 0s

    // SVD decomposition of F'
    svd_decompose(F_e, U_F, S_F, V_F);

//    std::cout << "S_F: " << S_F << std::endl;

    // enforce rank 2
    S_F(2, 2) = 0.0;
//    std::cout << "S_F rank 2: " << S_F << std::endl;

    Matrix<double> svd_res = U_F * S_F * V_F;
    mat3 F = to_mat3(svd_res);  // to fixed size for efficiency

    return F;
}


/// Step 2: Recover relative pose (R & t)
vec3 triangulate_linear(vec3 point0, vec3 point1, mat3 K, mat3 R, vec3 t){

//    std::cout << "R current:\n" << R << std::endl;
//    std::cout << "t current:\n" << t << std::endl;

    // P --> MP = | u |      P --> M'P = | u' |
    //            | v |                  | v' |
    //            | 1 |                  | 1  |

    // projection matrices of the two cameras (M & M' or M0 & M1 here)

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

//    vec3 point0_c = inverse(K) * point0;
//    vec3 point1_c = inverse(K) * point1;

//    std::cout << "point 0: " << point0 << std::endl;
//    std::cout << "point 1: " << point1 << std::endl;
//    std::cout << "point 0 canon: " << point0_c << std::endl;
//    std::cout << "point 1 canon: " << point1_c << std::endl;

    // fill A
    Matrix<double> A(4, 4, 0.0);
    vec4 elem00 = point0.x * M0.row(2) - M0.row(0);
    vec4 elem10 = point0.y * M0.row(2) - M0.row(1);
    vec4 elem20 = point1.x * M1.row(2) - M1.row(0);
    vec4 elem30 = point1.y * M1.row(2) - M1.row(1);

//    elem00 = elem00 / norm(elem00);
//    elem10 = elem10 / norm(elem10);
//    elem20 = elem20 / norm(elem20);
//    elem30 = elem30 / norm(elem30);

//    std::cout << "elem 00:\n" << elem00 << std::endl;
//    std::cout << "elem 10:\n" << elem10 << std::endl;
//    std::cout << "elem 20:\n" << elem20 << std::endl;
//    std::cout << "elem 30:\n" << elem30 << std::endl;

    for (int i = 0; i < 4; ++i) {
        A(0,i) = elem00[i];
        A(1,i) = elem10[i];
        A(2,i) = elem20[i];
        A(3,i) = elem30[i];
    }

//    std::cout << "A:\n" << A << std::endl;

    int h = A.rows();
    int w = A.cols();

    Matrix<double> U(h, h, 0.0);
    Matrix<double> S(h, w, 0.0);
    Matrix<double> V(w, w, 0.0);

    svd_decompose(A, U, S, V);

//    std::cout << "U (A):\n" << U << std::endl;
//    std::cout << "S (A):\n" << S << std::endl;
//    std::cout << "V (A):\n" << V << std::endl;

    Matrix<double> P(4, 1, 0.0);
    for (int i = 0; i < 4; ++i) {
        P(i,0) = V(i, w-1);
    }
//    std::cout << "P:\n" << P << std::endl;
//    std::cout << "P norm:\n" << P / P(3, 0) << std::endl;
//
//    std::cout << "A * P:\n" << A * P << std::endl;
//    std::cout << "A * P_norm:\n" << A * (P / P(3,0)) << std::endl;

    Matrix<double> res = (P / P(3,0));
    vec3 point_3d = {float(res(0,0)), float(res(1,0)), float(res(2,0))};

//    point_3d = K * point_3d;

//    vec4 point_3d_hom = {float(P(0,0)), float(P(1,0)), float(P(2,0)), float(P(3,0))};
//
//    vec3 backpt0 = M0 * point_3d_hom;
//    vec3 backpt1 = M1 * point_3d_hom;
//
//    std::cout << "point 0 current:\n" << point0 << std::endl;
//    std::cout << "point 1 current:\n" << point1 << std::endl;
//
//    std::cout << "M0 * P: " << backpt0.x / backpt0.z << " , " << backpt0.y / backpt0.z << std::endl;
//    std::cout << "M1 * P: " << backpt1.x / backpt1.z << " , " << backpt1.y / backpt1.z << std::endl;

    return point_3d;
}

vec3 triangulate_linear_print(vec3 point0, vec3 point1, mat3 K, mat3 R, vec3 t){

//    std::cout << "R current:\n" << R << std::endl;
//    std::cout << "t current:\n" << t << std::endl;

    // P --> MP = | u |      P --> M'P = | u' |
    //            | v |                  | v' |
    //            | 1 |                  | 1  |

    // projection matrices of the two cameras (M & M' or M0 & M1 here)

    // M = K[I 0]
    mat34 M0(1.0f);
    M0(2,3) = 0.0;
    M0 = K * M0;

    // M' = K'[R t]
    mat34 M1(1.0f);
    for (int i = 0; i < R.num_rows(); ++i) {
        vec3 row = R.row(i);
        M1.set_row(i, vec4(row.x, row.y, row.z, 0.0));
    }
    M1.set_col(3, t);
    M1 = K * M1;

    vec3 point0_c = inverse(K) * point0;
    vec3 point1_c = inverse(K) * point1;

    std::cout << "point 0: " << point0 << std::endl;
    std::cout << "point 0 canon: " << point0_c << std::endl;
    std::cout << "point 1: " << point1 << std::endl;
    std::cout << "point 1 canon: " << point1_c << std::endl;

    // fill A
    Matrix<double> A(4, 4, 0.0);
    vec4 elem00 = point0_c.x * M0.row(2) - M0.row(0);
    vec4 elem10 = point0_c.y * M0.row(2) - M0.row(1);
    vec4 elem20 = point1_c.x * M1.row(2) - M1.row(0);
    vec4 elem30 = point1_c.y * M1.row(2) - M1.row(1);

    elem00 = elem00 / norm(elem00);
    elem10 = elem10 / norm(elem10);
    elem20 = elem20 / norm(elem20);
    elem30 = elem30 / norm(elem30);

//    std::cout << "elem 00:\n" << elem00 << std::endl;
//    std::cout << "elem 10:\n" << elem10 << std::endl;
//    std::cout << "elem 20:\n" << elem20 << std::endl;
//    std::cout << "elem 30:\n" << elem30 << std::endl;

    for (int i = 0; i < 4; ++i) {
        A(0,i) = elem00[i];
        A(1,i) = elem10[i];
        A(2,i) = elem20[i];
        A(3,i) = elem30[i];
    }

    std::cout << "A:\n" << A << std::endl;

    int h = A.rows();
    int w = A.cols();

    Matrix<double> U(h, h, 0.0);
    Matrix<double> S(h, w, 0.0);
    Matrix<double> V(w, w, 0.0);

    svd_decompose(A, U, S, V);

    std::cout << "U (A):\n" << U << std::endl;
    std::cout << "S (A):\n" << S << std::endl;
    std::cout << "V (A):\n" << V << std::endl;

    Matrix<double> P(4, 1, 0.0);
    for (int i = 0; i < 4; ++i) {
        P(i,0) = V(i, w-1);
    }
    std::cout << "P:\n" << P << std::endl;
    std::cout << "P norm:\n" << P / P(3, 0) << std::endl;

    std::cout << "A * P:\n" << A * P << std::endl;
    std::cout << "A * P_norm:\n" << A * (P / P(3,0)) << std::endl;

    Matrix<double> res = (P / P(3,0));
    vec3 point_3d = {float(res(0,0)), float(res(1,0)), float(res(2,0))};

//    vec4 point_3d_hom = {float(P(0,0)), float(P(1,0)), float(P(2,0)), float(P(3,0))};
//    vec4 point_3d_hom = {float(res(0,0)), float(res(1,0)), float(res(2,0)), 1.0};
//
//    vec3 backpt0 = M0 * point_3d_hom;
//    vec3 backpt1 = M1 * point_3d_hom;
//
//    std::cout << "point 0 current:\n" << point0 << std::endl;
//    std::cout << "point 1 current:\n" << point1 << std::endl;
//
//    std::cout << "M0 * P: " << backpt0.x / backpt0.z << " , " << backpt0.y / backpt0.z << std::endl;
//    std::cout << "M1 * P: " << backpt1.x / backpt1.z << " , " << backpt1.y / backpt1.z << std::endl;

    Matrix<double> point_3d_hom(4,1);
    point_3d_hom.set_column({point_3d.x, point_3d.y, point_3d.z, 1.0}, 0);

    std::cout << "point 0 current:\n" << point0 << std::endl;
    std::cout << "point 1 current:\n" << point1 << std::endl;

    auto backpt0 = to_Matrix(M0) * point_3d_hom;
    auto backpt1 = to_Matrix(M1) * point_3d_hom;

    std::cout << "tested pt0: " << backpt0 << std::endl;
    std::cout << "tested pt1: " << backpt1 << std::endl;
    std::cout << "M0 * P: " << backpt0(0,0) / backpt0(2,0) << " , "
                            << backpt0(1,0) / backpt0(2,0) << std::endl;
    std::cout << "M1 * P: " << backpt1(0,0) / backpt1(2,0) << " , "
                            << backpt1(1,0) / backpt1(2,0) << std::endl;



    return point_3d;
}


std::tuple<mat3, vec3> correct_direction_poses(mat3 R01, mat3 R02, vec3 u3, mat3 K,
                                               std::vector<vec3> pts0,
                                               std::vector<vec3> pts1){
//    std::cout << "R1:\n" << R01 << std::endl;
//    std::cout << "R2:\n" << R02 << std::endl;
//
//    std::cout << "det R1:\n" << determinant(R01) << std::endl;
//    std::cout << "det R2:\n" << determinant(R02) << std::endl;
//
//    std::cout << "u3: " << u3 << std::endl;

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

//            std::cout << "resulting 3D point:\n" << pt_3d << std::endl;

            // find if the point is in front of both cameras
            // todo: camera center as simply t_curr.z was not correct
            vec3 camera_center_1 = R_curr * t_curr;
            if (pt_3d.z > 0 && pt_3d.z > camera_center_1.z){
                cnt ++;
            }
        }
        // if more points are in front of both cameras in this situation
        if (cnt > cnt_max){
            cnt_max = cnt;
            R_best = R_curr;
            t_best = t_curr;
        }
    }
    std::cout << "largest nr of valid points: " << cnt_max << std::endl;

//    std::cout << "det R1 corrected:\n" << determinant(R1) << std::endl;
//    std::cout << "det R2 corrected:\n" << determinant(R2) << std::endl;
//    std::cout << "det R best:\n" << determinant(R_best) << std::endl;

    std::cout << "t:\n" << t_best << std::endl;
    std::cout << "R:\n" << R_best << std::endl;

    return std::make_tuple(R_best, t_best);    // std::make_tuple(R2, t2);
};


std::tuple<mat3, vec3> relative_position(mat3 E, mat3 K, std::vector<vec3> points_0, std::vector<vec3> points_1){
    // determine relative position image 1 to image 0
    // as translation t and rotation R

    // W =  |   0  -1   0   |      Z =  |   0   1   0   |
    //      |   1   0   0   |           |  -1   0   0   |
    //      |   0   0   1   |           |   0   0   0   |

    mat3 W = {0, -1, 0, 1, 0, 0, 0, 0, 1};
    mat3 Z = {0, 1, 0, -1, 0, 0, 0, 0, 0};

    // E = U S V^T
    // E = [t]x R

    // [t]x = U Z U^T
    // t = +- u3
    // R = U W V^T or U W^T V^T

//    std::vector<vec3> points_0_canon;
//    std::vector<vec3> points_1_canon;
//
//    for (vec3 point : points_0){
//        points_0_canon.push_back(inverse(K) * point);
//    }
//
//    for (vec3 point : points_1){
//        points_1_canon.push_back(inverse(K) * point);
//    }


    // SVD decomposition of E
    Matrix<double> U_mat(3, 3, 0.0);
    Matrix<double> S_mat(3, 3, 0.0);
    Matrix<double> V_mat(3, 3, 0.0);

    Matrix<double> E_mat = to_Matrix(E);
    svd_decompose(E_mat, U_mat, S_mat, V_mat);

//    std::cout << "U (E):\n" << U_mat << std::endl;
//    std::cout << "S (E):\n" << S_mat << std::endl;
//    std::cout << "V^T (E):\n" << V_mat << std::endl;

//    std::cout << "E:\n" << E << std::endl;

//    S_mat(0,0) = 1.0;
//    S_mat(1,1) = 1.0;
//
//    std::cout << "S (E) ref:\n" << S_mat << std::endl;

//    Matrix<double> E_assume = U_mat * S_mat * V_mat;
//
//    Matrix<double> U_ass(3, 3, 0.0);
//    Matrix<double> S_ass(3, 3, 0.0);
//    Matrix<double> V_ass(3, 3, 0.0);
//
//    svd_decompose(E_assume, U_ass, S_ass, V_ass);

    mat3 U = to_mat3(U_mat);
    mat3 S = to_mat3(S_mat);
    mat3 V = to_mat3(V_mat);

    // rotation R
    mat3 R1 = U * W * V;
    mat3 R2 = U * transpose(W) * V;

    // translation t
    vec3 u3 = U.col(2);

    // R and t, both in the correct direction
    std::tuple<mat3, vec3> rel_pos_corr = correct_direction_poses(R1, R2, u3, K, points_0, points_1);

    return rel_pos_corr;
};


/// Step 3: Determine the 3D coordinates
void coords_3d(){};



/**
 * TODO: Finish this function for reconstructing 3D geometry from corresponding image points.
 * @return True on success, otherwise false. On success, the reconstructed 3D points must be written to 'points_3d'.
 */
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

    // TODO: check if the input is valid (always good because you never known how others will call your function).
    test_input();

    // TODO: Estimate relative pose of two views. This can be subdivided into:

    /// estimate the fundamental matrix F
    // normalize input points
    std::tuple<std::vector<vec3>, mat3 > norm0 = normalize_input(points_0);
    std::tuple<std::vector<vec3>, mat3 > norm1 = normalize_input(points_1);

    // get normalized points
    std::vector<vec3> pts0_norm = std::get<0>(norm0);
    std::vector<vec3> pts1_norm = std::get<0>(norm1);
    // get normalization transformation matrices
    mat3 T_0 = std::get<1>(norm0);
    mat3 T_1 = std::get<1>(norm1);

    // construct F
    mat3 F_init = estimate_F(pts0_norm, pts1_norm);
    // de-normalize F
    mat3 F_norm = transpose(T_1) * F_init * T_0;
    // rescale F
    mat3 F = F_norm / F_norm(2,2);
//    mat3 F = F_norm;
//    F(2,2) = 1.0;

//    mat3 F_nonorm = estimate_F(points_0, points_1);
//    F_nonorm = F_nonorm / F_nonorm(2,2);
//
//    F = F_nonorm;


    std::cout << "F:\n" << F << std::endl;

    Matrix<double> F_test = to_Matrix(F);
    Matrix<double> ptest0(1, 3, 0.0);
    ptest0.set_row({points_0[0].x, points_0[0].y, points_0[0].z}, 0);
    Matrix<double> ptest1(3, 1, 0.0);
    ptest1.set_column({points_1[0].x, points_1[0].y, points_1[0].z}, 0);

//    std::cout << "ptest0: " << ptest0 << std::endl;
//    std::cout << "ptest1: " << ptest1 << std::endl;

    Matrix<double> testF = ptest0 * F_test * ptest1 ;

//    std::cout << "test F:" << testF << std::endl;
//    std::cout << "test determinant F:" << determinant(F) << std::endl;


    /// compute the essential matrix E;
    // construct K from input intrinsic parameters
    // K =  |   fx    skew    cx   |
    //      |   0      fy     cy   |
    //      |   0      0       1   |
    mat3 K(1.0f);
    K(0,0) = fx;
    K(1,1) = fy;
    K(0,2) = cx;
    K(1,2) = cy;

    std::cout << "K:\n" << K << std::endl;

    // E = K^T F K  --> intrinsics are known
    // difference E and F: F is image coordinates, E = camera coordinates
    // K' = K according to assignment
    mat3 E = transpose(K) * F * K;

//    vec3 point0_c = inverse(K) * points_0[0];
//    vec3 point1_c = inverse(K) * points_1[0];
//
//    Matrix<double> ptest0e(1, 3, 0.0);
//    ptest0.set_row({point0_c.x, point0_c.y, point0_c.z}, 0);
//    Matrix<double> ptest1e(3, 1, 0.0);
//    ptest1.set_column({point1_c.x, point1_c.y, point1_c.z}, 0);
//
//    std::cout << "test E 2.0: " << ptest0e * E * ptest1e << std::endl;

    std::cout << "E:\n" << E << std::endl;

    std::cout << "test E:" <<  ptest0 * to_Matrix(E) * ptest1 << std::endl;
    std::cout << "test det E:" <<  determinant(E) << std::endl;
    std::cout << "test F with E:" <<  transpose(inverse(K)) * E * inverse(K) << std::endl;

    /// recover rotation R and t.
    std::tuple<mat3, vec3> pos = relative_position(E, K, points_0, points_1);
    R = std::get<0>(pos);
    t = std::get<1>(pos);

//    std::cout << "rr" << R * transpose(R) << std::endl;


    // TODO: Reconstruct 3D points. The main task is
    //      - triangulate a pair of image points (i.e., compute the 3D coordinates for each corresponding point pair)
    for (int i = 0; i < points_0.size(); ++i) {
        vec3 pt3d = triangulate_linear(points_0[i], points_1[i], K, R, t);

//        std::cout << "resulting 3d point: " << pt3d << std::endl;

        points_3d.push_back(pt3d);
    }


    // TODO: Don't forget to
    //          - write your recovered 3D points into 'points_3d' (the viewer can visualize the 3D points for you);
    //          - write the recovered relative pose into R and t (the view will be updated as seen from the 2nd camera,
    //            which can help you to check if R and t are correct).
    //       You must return either 'true' or 'false' to indicate whether the triangulation was successful (so the
    //       viewer will be notified to visualize the 3D points and update the view).
    //       However, there are a few cases you should return 'false' instead, for example:
    //          - function not implemented yet;
    //          - input not valid (e.g., not enough points, point numbers don't match);
    //          - encountered failure in any step.



    return points_3d.size() > 0;
}