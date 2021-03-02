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
 * TODO: Finish this function for calibrating a camera from the corresponding 3D-2D point pairs.
 *       You may define a few functions for some sub-tasks.
 *
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
    // todo: remove this in final version
    std::cout << std::endl;
    std::cout << "TODO: I am going to implement the calibration() function in the following file:" << std::endl
              << "\t" << __FILE__ << std::endl;
    std::cout << "TODO: After implementing the calibration() function, I will disable all unrelated output ...\n\n";

    /// TO DO: check if input is valid (e.g., number of correspondences >= 6, sizes of 2D/3D points must match)
    // Already partially implemented:
    //      - in open() method: when reading the file, if a line does not contain exactly 5 values,
    //        none of the coordinates are added to points_2d and points_3d.
    //      - in key_press_event() method: throws error if the size of points_2d or points_3d < 6.

    //check for duplicate and negative points in input:
    std::vector<int> duplicateLocations;

    for (int i = 0; i < points_3d.size(); i ++ ){
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
        for (int j = 0; j < points_3d.size(); j ++) {
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
    int height = points_2d.size() * 2;  // 2n
    Matrix<double> P(height, 12,  0.0);
    std::cout << "P: \n" << P << std::endl;

    // x_coor_pi * (m3 * Pi) - m1 * Pi = 0
    // y_coor_pi * (m3 * Pi) - m2 * Pi = 0

    // P matrix is entire system of equations (see above)
    // m = M as a vector of size (1, 12)
    // P size = (2n, (4 * 3))
    // P * m = 0

    // filling P: for each point pair
    for (int i = 0; i < points_2d_.size(); ++i) {
        // P_i = real world coordinate
        double Px = points_3d_[i][0];
        double Py = points_3d_[i][1];
        double Pz = points_3d_[i][2];
        double Pw = 1.0;    // homogenous coordinates

        // filling the rows of P: twice, for x and y of the 2d point (u & v)
        for (int j = 0; j < 2; ++j) {
            // Pi^T
            P(i*2 + j, 0 + j*4) = Px;
            P(i*2 + j, 1 + j*4) = Py;
            P(i*2 + j, 2 + j*4) = Pz;
            P(i*2 + j, 3 + j*4) = Pw;
            // 0^T
            P(i*2 + j, 4 - j*4) = 0.0;
            P(i*2 + j, 5 - j*4) = 0.0;
            P(i*2 + j, 6 - j*4) = 0.0;
            P(i*2 + j, 7 - j*4) = 0.0;
            // -[u or v] * Pi^T
            P(i*2 + j, 8) = -points_2d_[i][j] * Px;
            P(i*2 + j, 9) = -points_2d_[i][j] * Py;
            P(i*2 + j, 10) = -points_2d_[i][j] * Pz;
            P(i*2 + j, 11) = -points_2d_[i][j] * Pw;
        }
    }

    std::cout << "P: \n" << P << std::endl;
    // todo: does the SVD work with size (2n, 12) or will it need (2n, 4)?

    /// TASK: solve for M (the whole projection matrix, i.e., M = K * [R, t]) using SVD decomposition.
    ///   Optional: you can check if your M is correct by applying M on the 3D points. If correct, the projected point
    ///             should be very close to your input images points.

    Matrix<double> U(height, height, 0.0);   // initialized with 0s
    Matrix<double> S(height, 12, 0.0);   // initialized with 0s
    Matrix<double> V(12, 12, 0.0);   // initialized with 0s

    // Compute the SVD decomposition of P
    svd_decompose(P, U, S, V);

    // M from m (m = last column of V)
    Matrix<double> M(3, 4, 0.0);    // initialized with 0s

    for (int col = 0; col < 4; ++col) {
        for (int row = 0; row < 3; ++row) {
            M(row, col) = V(row + row*3 + col, 11);
        }
    }

    std::cout << "M: \n" << M << std::endl;

    // check if M is correct by applying it to the 3D points
    for (int i=0; i<points_2d_.size(); ++i) {
        std::vector<double> pts_3d = {points_3d_[i][0], points_3d_[i][1], points_3d_[i][2], 1.0};   // homogenous
        auto test_pts = M * pts_3d;
        std::cout << "\t real points: " << i << ": (" << points_3d_[i] << ") <-> (" << points_2d_[i] << ")" << std::endl;
        std::cout << "\t own points:  " << i << ": (" << round(test_pts[2]) << " "
                                                      << round(test_pts[1]) << " "
                                                      << round(test_pts[0]) << ")" << std::endl;
    }

    /// TASK: extract intrinsic parameters from M.

    /// TASK: extract extrinsic parameters from M.

    /// rotation matrix R

    // r1 = row 1 of R, etc.

    // r1 = (a2 x a3) / length(a2 x a3)



    /// translation 3D vector t (camera location)


    /// TASK: uncomment the line below to return true when testing your algorithm and in you final submission.
    // this draws a camera with the calculated M parameters
    return false;



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

















