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
    std::cout << std::endl; // extra space for clear output


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

    //construct the P matrix (so P * m = 0).

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

    // extract intrinsic parameters from M.
    // A = the three leftmost columns of M, b = column 4
    auto a1 = M.get_row(0);
    a1 = {a1[0], a1[1], a1[2]};
    auto a2 = M.get_row(1);
    a2 = {a2[0], a2[1], a2[2]};
    auto a3 = M.get_row(2);
    a3 = {a3[0], a3[1], a3[2]};

    // scaling factor rho
    double rho =  (1 / norm(a3) );
    // sign of rho
    if (M(2, 3) < 0){
        rho = rho * -1;
    }

    // principal point (cx, cy)
    auto u0 = pow(rho, 2) * a1 * a3;
    auto v0 = pow(rho, 2) * a2 * a3;

    cx = u0;
    cy = v0;

    // skew angle theta
    double theta = acos(- ( (cross(a1, a3) * cross(a2, a3)) / (norm(cross(a1, a3)) * norm(cross(a2, a3))) ));

    // focal length fx & fy
    double alpha = pow(rho, 2) * norm(cross(a1, a3)) * sin(theta);
    double beta = pow(rho, 2) * norm(cross(a2, a3)) * sin(theta);

    fx = (float) alpha;
    fy = (float) (beta / sin(theta));


    // skew factor
    skew = (float) (- fx * cos(theta));


    // rotation matrix R
    auto r1 = cross(a2, a3) / norm(cross(a2, a3));
    auto r3 = rho * a3;
    auto r2 = cross(r3, r1);

    R(0, 0) = r1[0];
    R(0, 1) = r1[1];
    R(0, 2) = r1[2];

    R(1, 0) = r2[0];
    R(1, 1) = r2[1];
    R(1, 2) = r2[2];

    R(2, 0) = r3[0];
    R(2, 1) = r3[1];
    R(2, 2) = r3[2];

    // translation 3D vector t (camera location)
    Matrix<double> K(3, 3, 0.0);   // initialized with 0s

    K(0, 0) = fx;
    K(0, 1) = skew;
    K(0, 2) = cx;
    K(1, 1) = fy;
    K(1, 2) = cy;
    K(2, 2) = 1.0;


    Matrix<double> invK(3, 3);
    inverse(K, invK);

    auto b_T = M.get_column(3);
    Matrix<double> b(3, 1, 0.0);
    b[0][0] = b_T[0];
    b[1][0] = b_T[1];
    b[2][0] = b_T[2];

    Matrix<double> t_T = rho * invK * b;

    t[0] = t_T[0][0];
    t[1] = t_T[1][0];
    t[2] = t_T[2][0];

    return true;
}

















