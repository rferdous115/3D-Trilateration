#include <iostream>
#include <iomanip>
#include <fstream>

#include "Eigen/Dense"

using namespace Eigen;
using namespace std;

int main() {
    double Ax = 103.;
    double Ay = 350.;
    double Az = 142.;

    double Bx = -127.;
    double By = 350.;
    double Bz = 128.;

    double Cx = -15.;
    double Cy = 500.;
    double Cz = 540.;

    double rDA = 172.064;
    double rDB = 183.527;
    double rDC = 401.103;

    // Method 4

    // Aim: to solve for Dx, Dy, and Dz, which are components of the vectors, and spheres:
    // |R_DA|^2 = (Dx - Ax)^2 + (Dy - Ay)^2 + (Dz - Az)^2
    // |R_DB|^2 = (Dx - Bx)^2 + (Dy - By)^2 + (Dz - Bz)^2
    // |R_DC|^2 = (Dx - Cx)^2 + (Dy - Cy)^2 + (Dz - Cz)^2

    // Preliminary Calculations

    // Let P1 be the vector from the center of sphere R_DA
    // Let P2 be the vector from the center of sphere R_DB
    // Let P3 be the vector from the center of sphere R_DC

    Vector3d P1(Ax, Ay, Az);
    Vector3d P2(Bx, By, Bz);
    Vector3d P3(Cx, Cy, Cz);


    cout << "P1" << endl;
    cout << P1 << endl;
    cout << endl;

    cout << "P2" << endl;
    cout << P2 << endl;
    cout << endl;

    cout << "P3" << endl;
    cout << P3 << endl;
    cout << endl;

    // e_x = (P1 - P2) / ||P2 - P1|| is the unit vector in the direction from P1 to P2.

    Vector3d temp1 = P2 - P1;
    Vector3d e_x = (P2 - P1) / temp1.norm();
    cout << "e_x" << endl;
    cout << e_x << endl;
    cout << endl;

    // i = e_x $dot$ (P3 - P1) is the signed magnitude of the x component, in the figure 1 component system, of the vector from P1 to P3.

    Vector3d temp2 = P3 - P1;

    float i = e_x.dot(temp2);
    cout << "i" << endl;
    cout << i << endl;
    cout << endl;

    // e_y = (P3 - P1 - (i * e_x)) / ||P3 - P1 - (i * e_x)|| is the unit vector in the y direction.
    // Note that the points P1, P2, and P3 are all in the z = 0 plane of the figure 1 coordinate system.

    Vector3d temp3 = temp2 - i * e_x;

    Vector3d e_y = (temp3)/temp3.norm();
    cout << "e_y" << endl;
    cout << e_y << endl;
    cout << endl;

    // e_z = e_x $cross$ e_y is the third basis unit vector. Therefore,

    Vector3d e_z = e_x.cross(e_y);
    cout << "e_z" << endl;
    cout << e_z << endl;
    cout << endl;

    // d = ||P2 - P1|| is the distance between the centers P1 and P2 ...

    float d = (P2 - P1).norm();
    cout << "d" << endl;
    cout << d << endl;
    cout << endl;

    // and j = e_y $dot$ (P3 - P1) is the signed magnitude of the y component, in the figure 1 coordinate system, of the vector from P1 to P3.

    float j = e_y.dot(temp2);
    cout << "j" << endl;
    cout << j << endl;
    cout << endl;

    // using i, d, and j, we will solve for x, y, and z in the Derivation section.

    // Derivation

    // For simplification of the calculations:
    // 1. all three centers are in the plane z = 0.
    // 2. the sphere center, P1, is at the origin, and
    // 3. the sphere center, P2, is on the x-axis.

    // We start with the equations for the three spheres:
    // r^2_1 = x^2 + y^2 + z^2
    // r^2_2 = (x - d)^2 + y^2 + z^2
    // r^2_3 = (x - i)^2 + (y - j)^2 + z^2

    float x;
    float y;
    float z;

    // solving for x:

    x = (rDA * rDA - rDB * rDB + d*d) / (2 * d);

    // solving for y:

    y = (rDA * rDA - rDC * rDC - 2 * i * x + i * i + j * j) / (2 * j);

    // Now that the x- and y-coordinates of the solution point are found, the formula can be rearranged for the first sphere to find the z-coordinate:

    float temp4 = rDA * rDA - x * x - y * y;

    z = sqrt(temp4);

    // Final Calculations:

    // Using i, d, and j as computed above, solve for x, y and z as described in the Derivation section. Then ...
    // p_1_2 = P1 + (x * e_x) + (y * e_y) +- (z * e_z)

    Vector3d p_1_2_sol1 = P1 + (x * e_x) + (y * e_y) + (z * e_z);
    Vector3d p_1_2_sol2 = P1 + (x * e_x) + (y * e_y) - (z * e_z);

    // Print solutions:

    cout << p_1_2_sol1 << endl;
    cout << endl;
    cout << p_1_2_sol2 << endl;

    // Success
    // Output:

    //    -3.40289
    //    485.194
    //    139.338
    //
    //    -8.7031
    //    249.986
    //    226.413

    return 0;
}
