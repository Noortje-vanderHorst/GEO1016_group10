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

#include "camera_calibration.h"
#include "matrix_algo.h"


using namespace easy3d;


/**
 * @param points_3d   An array of 3D points.
 * @param points_2d   An array of 2D points.
 * @return True on success, otherwise false. On success, the camera parameters are returned by
 *           - fx and fy: the focal length (in our slides, we use 'alpha' and 'beta'),
 *           - cx and cy: the principal point (in our slides, we use 'u0' and 'v0'),
 *           - skew:      the skew factor ('-alpha * cot_theta')
 *           - R:         the 3x3 rotation matrix encoding camera orientation.
 *           - t:         a 3D vector encoding camera location.
 */

bool CameraCalibration::calibration(
        const std::vector<vec3>& points_3d,
        const std::vector<vec2>& points_2d,
        float& fx, float& fy,
        float& cx, float& cy,
        float& skew,
        mat3& R,
        vec3& t)
{
    //check for duplicate and negative points in input:
    std::vector<int> duplicateLocations;

    for (unsigned int i = 0; i < points_3d.size(); i ++ ){
        // check for negative 3D points
        if (points_3d[i].x < 0 || points_3d[i].y < 0 || points_3d[i].x < 0){
            std::cout << "Invalid 3d point with negative coordinates: ("
                      << points_3d[i].x << " "
                      << points_3d[i].y << " "
                      << points_3d[i].z<< ")" << std::endl;
            std::cout << "Point is ignored" << std::endl;

            duplicateLocations.emplace_back(i);
            continue;
        }

        // check for negative 2D points
        if (points_2d[i].x < 0 || points_2d[i].y < 0){
            std::cout << "Invalid 2d point with negative coordinates: ("
                      << points_2d[i].x << " "
                      << points_2d[i].y << ")" << std::endl;
            std::cout << "Point is ignored" << std::endl;

            duplicateLocations.emplace_back(i);
            continue;
        }

        // check for duplicates
        for (unsigned int j = 0; j < points_3d.size(); j ++) {
            if ( i >= j) {
                continue;
            }
            if (points_3d[i].x == points_3d[j].x
                && points_3d[i].y == points_3d[j].y
                && points_3d[i].z == points_3d[j].z) {

                std::cout << "Duplicate 3d point with coordinates: ("
                          << points_3d[i][0] << " "
                          << points_3d[i][1] << " "
                          << points_3d[i][2] << ")" << std::endl;
                std::cout << "Point is ignored" << std::endl;

                duplicateLocations.emplace_back(i);
                break;
            }
        }
    }

    // check if the size after removal is big enough
    if (points_3d.size() - duplicateLocations.size() < 6){
        std::cerr << "expecting at least 6 unique pairs of 3D/2D corresponding points" << std::endl;
        return false;
    }

    // remove invalid points from data:
    for(int i = duplicateLocations.size(); i--;){
        points_3d_.erase(points_3d.begin() + i);
        points_2d_.erase(points_2d.begin() + i);
    }

    /// TASK: construct the P matrix (so P * m = 0).

    // P initialized with 0s, size (2*[number of point pairs], 3)
    int height = (int) points_2d.size() * 2;  // 2n
    Matrix<double> P(height, 12,  0.0);

    // filling P: for each point pair
    for (int i = 0; i < (int) points_2d_.size(); ++i) {
        // P_i = real world coordinate
        double Px = points_3d_[i][0];
        double Py = points_3d_[i][1];
        double Pz = points_3d_[i][2];

        // filling the rows of P: twice, for x and y of the 2d point (u & v)
        for (int j = 0; j < 2; ++j) {
            // Pi^T
            P(i*2 + j, 0 + j*4) = Px;
            P(i*2 + j, 1 + j*4) = Py;
            P(i*2 + j, 2 + j*4) = Pz;
            P(i*2 + j, 3 + j*4) = 1.0;
            // -[u or v] * Pi^T
            P(i*2 + j, 8) = -points_2d_[i][j] * Px;
            P(i*2 + j, 9) = -points_2d_[i][j] * Py;
            P(i*2 + j, 10) = -points_2d_[i][j] * Pz;
            P(i*2 + j, 11) = -points_2d_[i][j] * 1.0;
        }
    }

    std::cout << "P: \n" << P << std::endl;

    //solve for M (the whole projection matrix, i.e., M = K * [R, t]) using SVD decomposition.

    Matrix<double> U(height, height, 0.0);   // initialized with 0s
    Matrix<double> S(height, 12, 0.0);   // initialized with 0s
    Matrix<double> V(12, 12, 0.0);   // initialized with 0s

    // Compute the SVD decomposition of P
    svd_decompose(P, U, S, V);

    // M from m (m = last column of V)
    Matrix<double> M(3, 4, 0.0);    // initialized with 0s

    // populate M
    for (int i = 0; i < 3; i ++){
        for (int j = 0; j < 4; j++){
            M(i,j) = V(i * 4 + j, 11);
        }
    }
    std::cout << "M: \n" << M << std::endl;

    // check if M is correct by applying it to the 3D points
    for (int i=0; i<points_2d_.size(); ++i) {
        Matrix<double> pts_3d(4, 1, 0.0);
        pts_3d[0][0] = points_3d_[i][0];
        pts_3d[1][0] = points_3d_[i][1];
        pts_3d[2][0] = points_3d_[i][2];
        pts_3d[3][0] = 1.0;
        auto test_pts = M * pts_3d;
        std::cout << "\t real points: " << i << ": (" << points_3d_[i] << ") <-> (" << points_2d_[i] << ")" << std::endl;
        std::cout << "\t own points:  " << i << ": (" << test_pts[0][0] / test_pts[2][0] << " "
                                                      << test_pts[0][1] / test_pts[2][0] << ")" << std::endl;

    }

    // extract intrinsic parameters from M.
    // A = the three leftmost columns of M, b = column 4
    auto a1 = M.get_row(0);
    a1 = {a1[0], a1[1], a1[2]};
    auto a2 = M.get_row(1);
    a2 = {a2[0], a2[1], a2[2]};
    auto a3 = M.get_row(2);
    a3 = {a3[0], a3[1], a3[2]};

    std::cout << "A_1: " << a1 << std::endl;
    std::cout << "A_2: " << a2 << std::endl;
    std::cout << "A_3: " << a3 << std::endl;

    /// scaling factor rho
    double rho =  1 / norm(a3);

    std::cout << "rho: " << rho << std::endl;
    std::cout << "norm(a3): " << norm(a3) << std::endl;
    std::cout << "M: " << M << std::endl;

    /// principal point (cx, cy)
    auto u0 = pow(rho, 2) * a1 * a3;
    auto v0 = pow(rho, 2) * a2 * a3;

    std::cout << "cx: " << u0 << std::endl;
    std::cout << "cy: " << v0 << std::endl;
    cx = u0;
    cy = v0;

    /// skew angle theta
    double theta = acos(- ( (cross(a1, a3) * cross(a2, a3)) / (norm(cross(a1, a3)) * norm(cross(a2, a3))) ));

    std::cout << "theta: " << theta << std::endl;
    std::cout << "theta (deg): " << rad2deg(theta) << std::endl;

    /// focal length fx & fy

    // alpha = rho^2 * length(a1 x a3) * sin(theta)
    // beta = rho^2 * length(a2 x a3) * sin(theta)

    double alpha = pow(rho, 2) * norm(cross(a1, a3)) * sin(theta);
    double beta = pow(rho, 2) * norm(cross(a2, a3)) * sin(theta);
    std::cout << "alpha (fx): " << alpha << std::endl;
    std::cout << "beta: " << beta << std::endl;
    std::cout << "fy: " << beta / sin(theta) << std::endl;

    fx = (float) alpha;
    fy = (float) (beta / sin(theta));


    /// skew factor
    skew = (float) (- fx * cos(theta));
    std::cout << "skew factor: " << skew << std::endl;


    /// rotation matrix R

    auto r1 = cross(a2, a3) / norm(cross(a2, a3));
    auto r3 = rho * a3;
    auto r2 = cross(r3, r1);

    std::cout << "r1: " << r1 << std::endl;
    std::cout << "r2: " << r2 << std::endl;
    std::cout << "r3: " << r3 << std::endl;

//    std::cout << "R: " << R << std::endl;

    R(0, 0) = r1[0];
    R(0, 1) = r1[1];
    R(0, 2) = r1[2];

    R(1, 0) = r2[0];
    R(1, 1) = r2[1];
    R(1, 2) = r2[2];

    R(2, 0) = r3[0];
    R(2, 1) = r3[1];
    R(2, 2) = r3[2];

    std::cout << "R: " << R << std::endl;

    /// translation 3D vector t (camera location)

    // t = rho * K^-1 * b

    // K =  | fx    s     cx  |
    //      | 0     fy    cy  |
    //      | 0     0     1   |

    Matrix<double> K(3, 3, 0.0);   // initialized with 0s

    K(0, 0) = fx;
    K(0, 1) = skew;
    K(0, 2) = cx;
    K(1, 1) = fy;
    K(1, 2) = cy;
    K(2, 2) = 1.0;

    std::cout << "K: " << K << std::endl;

    Matrix<double> invK(3, 3);
    inverse(K, invK);

    // todo: is b then the last column of M?
    auto b_T = M.get_column(3);
    Matrix<double> b(3, 1, 0.0);
    b[0][0] = b_T[0];
    b[1][0] = b_T[1];
    b[2][0] = b_T[2];


    Matrix<double> t_T = rho * invK * b;

    std::cout << "inv K: " << invK << std::endl;
    std::cout << "K * inv K: " << K * invK << std::endl;
    std::cout << "t_T: " << t_T << std::endl;

    t[0] = t_T[0][0];
    t[1] = t_T[1][0];
    t[2] = t_T[2][0];

    std::cout << "t: " << t << std::endl;


    /// TASK: uncomment the line below to return true when testing your algorithm and in you final submission.
    // this draws a camera with the calculated M parameters
    return true;



    /// TASK: The following code is just an example showing you SVD decomposition, matrix inversion, and some related.
    /// TASK: Delete the code below (or change "#if 1" in the first line to "#if 0") in you final submission.
#if 0
    std::cout << "[Liangliang:] Camera calibration requires computing the SVD and inverse of matrices.\n"
                 "\tIn this assignment, I provide you with a Matrix data structure for storing matrices of arbitrary\n"
                 "\tsizes (see matrix.h). I also wrote the example code to show you how to:\n"
                 "\t\t- use the dynamic 1D array data structure 'std::vector' from the standard C++ library;\n"
                 "\t\t  The points (both 3D and 2D) are stored in such arrays;\n"
                 "\t\t- use the template matrix class (which can have an arbitrary size);\n"
                 "\t\t- compute the SVD of a matrix;\n"
                 "\t\t- compute the inverse of a matrix;\n"
                 "\t\t- compute the transpose of a matrix.\n"
                 "\tThe following are just the output of these examples. You should delete ALL unrelated code and\n"
                 "\tavoid unnecessary output in you final submission.\n\n";

    // This is a 1D array of 'double' values. Alternatively, you can use 'double mat[25]' but you cannot change it
    // length. With 'std::vector', you can do append/delete/insert elements, and much more. The 'std::vector' can store
    // not only 'double', but also any other types of objects. In case you may want to learn more about 'std::vector'
    // check here: https://en.cppreference.com/w/cpp/container/vector
    std::vector<double> array = {1, 3, 3, 4, 7, 6, 2, 8, 2, 8, 3, 2, 4, 9, 1, 7, 3, 23, 2, 3, 5, 2, 1, 5, 8, 9, 22};
    array.push_back(5); // append 5 to the array (so the size will increase by 1).
    array.insert(array.end(), 10, 3);  // append ten 3 (so the size will grow by 10).

    // To access its values
    for (int i=0; i<array.size(); ++i)
        std::cout << array[i] << " ";  // use 'array[i]' to access its i-th element.
    std::cout << std::endl;

    // Define an m-by-n double valued matrix.
    // Here I use the above array to initialize it. You can also use A(i, j) to initialize/modify/access its elements.
    const int m = 6, n = 5;
    Matrix<double> A(m, n, array.data());    // 'array.data()' returns a pointer to the array.
    std::cout << "M: \n" << A << std::endl;

    Matrix<double> U(m, m, 0.0);   // initialized with 0s
    Matrix<double> S(m, n, 0.0);   // initialized with 0s
    Matrix<double> V(n, n, 0.0);   // initialized with 0s

    // Compute the SVD decomposition of A
    svd_decompose(A, U, S, V);

    // Now let's check if the SVD result is correct

    // Check 1: U is orthogonal, so U * U^T must be identity
    std::cout << "U*U^T: \n" << U * transpose(U) << std::endl;

    // Check 2: V is orthogonal, so V * V^T must be identity
    std::cout << "V*V^T: \n" << V * transpose(V) << std::endl;

    // Check 3: S must be a diagonal matrix
    std::cout << "S: \n" << S << std::endl;

    // Check 4: according to the definition, A = U * S * V^T
    std::cout << "M - U * S * V^T: \n" << A - U * S * transpose(V) << std::endl;

    // Define a 5 by 5 square matrix and compute its inverse.
    Matrix<double> B(5, 5, array.data());    // Here I use part of the above array to initialize B
    // Compute its inverse
    Matrix<double> invB(5, 5);
    inverse(B, invB);
    // Let's check if the inverse is correct
    std::cout << "B * invB: \n" << B * invB << std::endl;

    return false;
    // TODO: delete the above code in you final submission (which are just examples).
#endif
}

















