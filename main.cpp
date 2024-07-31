//
//  main.cpp
//  ffd
//
//  Created by Gary Chu on 28/02/2019.
//  Copyright © 2019 Gary Chu. All rights reserved.
//
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <OpenGL/gl.h>
#include <OpenGl/glu.h>

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <iterator>
#include "Eigen/Dense"

using namespace std;
using namespace Eigen;



MatrixXd Vertex_coor(214,2);
MatrixXd Face(332,3);

MatrixXd x_coor(214,1);
MatrixXd y_coor(214,1);
//double x_coor[214], y_coor[214];
int f_1[332], f_2[332], f_3[332];
int v_count = 0, f_count = 0;
int mouse_click_mode = 0;
int control_point[214];
MatrixXd control_coor(1,2);
int num_control = 0;
int current_control = -1;
MatrixXd edge(545,4);
MatrixXd edge_vector(545,2);
MatrixXd GGG_value[545];
MatrixXd H_value[545];
int w = 1000;
MatrixXd T_value[545];

MatrixXd v_first;


int find_close_point(double x, double y) {
    int index = 0;
    double min = 100000.0;
    for (int i = 0; i < 214; i++) {
        double diff = pow((x_coor(i) - x), 2.0) + pow((y_coor(i) - y), 2.0);
        if (diff < min) {
            index = i;
            min = diff;
        }
    }
    return index;
}

void scale_adj() {
    MatrixXd A_2(edge.rows()+num_control,214);
    MatrixXd b_2_x(edge.rows()+num_control,1);
    MatrixXd b_2_y(edge.rows()+num_control,1);
    A_2 << MatrixXd::Zero(A_2.rows(),A_2.cols());
    for (int r = 0; r < edge.rows(); r++) {
        MatrixXd e_temp(2,1);
        e_temp << edge_vector(r,0),
                  edge_vector(r,1);
        MatrixXd b_temp = T_value[r] * e_temp;
        b_2_x(r) = b_temp(0);
        b_2_y(r) = b_temp(1);
        //j ==> 1
        //i ==> -1
        A_2(r,edge(r,1)) = 1;
        A_2(r,edge(r,0)) = -1;
    }
    int ab_pos = 0;
    for (int i = 0; i < 214; i++) {
        if (control_point[i] == 1) {
//            b_2_x(edge.rows()+ab_pos) = w * x_coor(i);
//            b_2_y(edge.rows()+ab_pos) = w * y_coor(i);
            b_2_x(edge.rows()+ab_pos) = w * v_first(2*i);
            b_2_y(edge.rows()+ab_pos) = w * v_first(2*i+1);
            A_2(edge.rows()+ab_pos,i) = w;
            ab_pos++;
        }
    }
    MatrixXd AAA = (A_2.transpose() * A_2).inverse() * A_2.transpose();
    MatrixXd v_sec_x = AAA * b_2_x;
    MatrixXd v_sec_y = AAA * b_2_y;
    for (int i = 0; i < 214; i++) {
        x_coor(i) = v_sec_x(i);
        y_coor(i) = v_sec_y(i);
    }
//    x_coor = AAA * b_2_x;
//    y_coor = AAA * b_2_y;
}

void calculate_T() {
    for (int r = 0; r < edge.rows(); r++) {
        if (edge(r,3) == -1) {
            MatrixXd v_temp(6,1);
//            v_temp << x_coor(int(edge(r,0))),
//                      y_coor(int(edge(r,0))),
//                      x_coor(int(edge(r,1))),
//                      y_coor(int(edge(r,1))),
//                      x_coor(int(edge(r,2))),
//                      y_coor(int(edge(r,2)));
            v_temp << v_first(2*int(edge(r,0))),
            v_first(2*int(edge(r,0))+1),
            v_first(2*int(edge(r,1))),
            v_first(2*int(edge(r,1))+1),
            v_first(2*int(edge(r,2))),
            v_first(2*int(edge(r,2))+1);
            MatrixXd g_temp(2,6);
            g_temp << GGG_value[r](0,0), GGG_value[r](0,1), GGG_value[r](0,2), GGG_value[r](0,3), GGG_value[r](0,4), GGG_value[r](0,5),
                      GGG_value[r](1,0), GGG_value[r](1,1), GGG_value[r](1,2), GGG_value[r](1,3), GGG_value[r](1,4), GGG_value[r](1,5);
            MatrixXd cs = g_temp * v_temp;
            MatrixXd t_temp(2,2);
            t_temp << cs(0), cs(1),
                     -cs(1), cs(0);
            t_temp *= 1 / (sqrt(pow(cs(0),2.0) + pow(cs(1),2.0)));
            T_value[r] = t_temp;
        } else {
            MatrixXd v_temp(8,1);
//            v_temp << x_coor(int(edge(r,0))),
//                      y_coor(int(edge(r,0))),
//                      x_coor(int(edge(r,1))),
//                      y_coor(int(edge(r,1))),
//                      x_coor(int(edge(r,2))),
//                      y_coor(int(edge(r,2))),
//                      x_coor(int(edge(r,3))),
//                      y_coor(int(edge(r,3)));
            v_temp << v_first(2*int(edge(r,0))),
            v_first(2*int(edge(r,0))+1),
            v_first(2*int(edge(r,1))),
            v_first(2*int(edge(r,1))+1),
            v_first(2*int(edge(r,2))),
            v_first(2*int(edge(r,2))+1),
            v_first(2*int(edge(r,3))),
            v_first(2*int(edge(r,3))+1);
            MatrixXd cs = GGG_value[r] * v_temp;
            MatrixXd t_temp(2,2);
            t_temp << cs(0), cs(1),
                     -cs(1), cs(0);
            t_temp *= 1 / (sqrt(pow(cs(0),2.0) + pow(cs(1),2.0)));
            T_value[r] = t_temp;
        }
    }
}

void sim_trans() {
    MatrixXd A_1(2*edge.rows()+2*num_control,2*214);
    A_1 << MatrixXd::Zero(A_1.rows(),A_1.cols());
    MatrixXd b_1(2*edge.rows()+2*num_control,1);
    b_1 << MatrixXd::Zero(b_1.rows(),b_1.cols());
    for (int r = 0; r < edge.rows(); r++) {
        //h00(2r,2a),h01(2r,2a+1),h10(2r+1,2a),h11(2r+1,2a+1)
        int a = edge(r,0);
        //h02(2r,2b),h03(2r,2b+1),h12(2r+1,2b),h13(2r+1,2b+1)
        int b = edge(r,1);
        //h04(2r,2c),h05(2r,2c+1),h14(2r+1,2c),h15(2r+1,2c+1)
        int c = edge(r,2);
        //h06(2r,2d),h07(2r,2d+1),h16(2r+1,2d),h17(2r+1,2d+1)
        int d = edge(r,3);
        
        A_1(2*r,2*a) = H_value[r](0,0);
        A_1(2*r,2*a+1) = H_value[r](0,1);
        A_1(2*r+1,2*a) = H_value[r](1,0);
        A_1(2*r+1,2*a+1) = H_value[r](1,1);
        
        A_1(2*r,2*b) = H_value[r](0,2);
        A_1(2*r,2*b+1) = H_value[r](0,3);
        A_1(2*r+1,2*b) = H_value[r](1,2);
        A_1(2*r+1,2*b+1) = H_value[r](1,3);
        
        A_1(2*r,2*c) = H_value[r](0,4);
        A_1(2*r,2*c+1) = H_value[r](0,5);
        A_1(2*r+1,2*c) = H_value[r](1,4);
        A_1(2*r+1,2*c+1) = H_value[r](1,5);
        
        if (d != -1){
            A_1(2*r,2*d) = H_value[r](0,6);
            A_1(2*r,2*d+1) = H_value[r](0,7);
            A_1(2*r+1,2*d) = H_value[r](1,6);
            A_1(2*r+1,2*d+1) = H_value[r](1,7);
        }
    }
    int b_pos = 0;
    for (int i = 0; i < 214; i++) {
        if (i == current_control) {
            A_1(2*edge.rows()+b_pos,2*i) = w;
            b_1(2*edge.rows()+b_pos) = w * control_coor(0);
            b_pos++;
            A_1(2*edge.rows()+b_pos,2*i+1) = w;
            b_1(2*edge.rows()+b_pos) = w * control_coor(1);
            b_pos++;
        } else if (control_point[i] == 1) {
            A_1(2*edge.rows()+b_pos,2*i) = w;
            b_1(2*edge.rows()+b_pos) = w * x_coor(i);
            b_pos++;
            A_1(2*edge.rows()+b_pos,2*i+1) = w;
            b_1(2*edge.rows()+b_pos) = w * y_coor(i);
            b_pos++;
        }
    }
    v_first = (A_1.transpose() * A_1).inverse() * A_1.transpose() * b_1;
//    for (int i = 0; i < 214; i++) {
//        x_coor(i) = v_first(2*i);
//        y_coor(i) = v_first(2*i+1);
//    }
}

void calculate_H() {
    for (int i = 0; i < 545; i++) {
        MatrixXd H(2,8);
        MatrixXd E(2,2);
        H << -1, 0, 1, 0, 0, 0, 0, 0,
             0, -1, 0, 1, 0, 0, 0, 0;
        E << edge_vector(i,0), edge_vector(i,1),
             edge_vector(i,1), -edge_vector(i,0);
        H -= E * GGG_value[i];
        if (GGG_value[i](0,6) == -1 & GGG_value[i](0,7) == -1 & GGG_value[i](1,6) == -1 & GGG_value[i](1,7) == -1) {
            H(0,6) = 0;
            H(0,7) = 0;
            H(1,6) = 0;
            H(0,7) = 0;
        }
        H_value[i] = H;
    }
}

void calculate_G() {
    for (int i = 0; i < 545; i++) {
        if (edge(i,3) == -1) {
            MatrixXd G(6,4);
            G << x_coor(int(edge(i,0))), y_coor(int(edge(i,0))), 1, 0,
                 y_coor(int(edge(i,0))), -x_coor(int(edge(i,0))), 0, 1,
                 x_coor(int(edge(i,1))), y_coor(int(edge(i,1))), 1, 0,
                 y_coor(int(edge(i,1))), -x_coor(int(edge(i,1))), 0, 1,
                 x_coor(int(edge(i,2))), y_coor(int(edge(i,2))), 1, 0,
                 y_coor(int(edge(i,2))), -x_coor(int(edge(i,2))), 0, 1;
            MatrixXd temp = (G.transpose() * G).inverse() * G.transpose();
            MatrixXd fin(2,8);
            fin << temp(0,0), temp(0,1), temp(0,2), temp(0,3), temp(0,4), temp(0,5), -1, -1,
                   temp(1,0), temp(1,1), temp(1,2), temp(1,3), temp(1,4), temp(1,5), -1, -1;
            GGG_value[i] = fin;
        } else {
            MatrixXd G(8,4);
            G << x_coor(int(edge(i,0))), y_coor(int(edge(i,0))), 1, 0,
                 y_coor(int(edge(i,0))), -x_coor(int(edge(i,0))), 0, 1,
                 x_coor(int(edge(i,1))), y_coor(int(edge(i,1))), 1, 0,
                 y_coor(int(edge(i,1))), -x_coor(int(edge(i,1))), 0, 1,
                 x_coor(int(edge(i,2))), y_coor(int(edge(i,2))), 1, 0,
                 y_coor(int(edge(i,2))), -x_coor(int(edge(i,2))), 0, 1,
                 x_coor(int(edge(i,3))), y_coor(int(edge(i,3))), 1, 0,
                 y_coor(int(edge(i,3))), -x_coor(int(edge(i,3))), 0, 1;
            MatrixXd temp = (G.transpose() * G).inverse() * G.transpose();
            MatrixXd fin(2,8);
            fin << temp(0,0), temp(0,1), temp(0,2), temp(0,3), temp(0,4), temp(0,5), temp(0,6), temp(0,7),
                   temp(1,0), temp(1,1), temp(1,2), temp(1,3), temp(1,4), temp(1,5), temp(1,6), temp(1,7);
            GGG_value[i] = fin;
        }
    }
}

void find_edge() {
    for (int r = 0; r < edge.rows(); r++) {
        for (int c = 0; c < edge.cols(); c++) {
            edge(r,c) = -1;
        }
    }
    for (int i = 0; i < 332; i++) {
        int smaller = -1;
        int greater = -1;
        if (f_1[i] > f_2[i]) {
            smaller = f_2[i];
            greater = f_1[i];
        } else {
            smaller = f_1[i];
            greater = f_2[i];
        }
        for (int r = 0; r < edge.rows(); r++) {
            if (edge(r,0) == -1) {
                edge(r,0) = smaller;
                edge(r,1) = greater;
                edge(r,2) = f_3[i];
                break;
            } else if (edge(r,0) == smaller) {
                if (edge(r,1) == greater) {
                    edge(r,3) = f_3[i];
                    break;
                } else {
                    continue;
                }
            } else {
                continue;
            }
        }
        if (f_1[i] > f_3[i]) {
            smaller = f_3[i];
            greater = f_1[i];
        } else {
            smaller = f_1[i];
            greater = f_3[i];
        }
        for (int r = 0; r < edge.rows(); r++) {
            if (edge(r,0) == -1) {
                edge(r,0) = smaller;
                edge(r,1) = greater;
                edge(r,2) = f_2[i];
                break;
            } else if (edge(r,0) == smaller) {
                if (edge(r,1) == greater) {
                    edge(r,3) = f_2[i];
                    break;
                } else {
                    continue;
                }
            } else {
                continue;
            }
        }
        if (f_3[i] > f_2[i]) {
            smaller = f_2[i];
            greater = f_3[i];
        } else {
            smaller = f_3[i];
            greater = f_2[i];
        }
        for (int r = 0; r < edge.rows(); r++) {
            if (edge(r,0) == -1) {
                edge(r,0) = smaller;
                edge(r,1) = greater;
                edge(r,2) = f_1[i];
                break;
            } else if (edge(r,0) == smaller) {
                if (edge(r,1) == greater) {
                    edge(r,3) = f_1[i];
                    break;
                } else {
                    continue;
                }
            } else {
                continue;
            }
        }
    }
//    int index = 0;
//    for (int r  = 0; r < edge.rows(); r++) {
//        if (edge(r,0) == -1) {
//            index = r;
//            break;
//        }
//    }
    for (int r = 0; r < edge.rows(); r++) {
        edge_vector(r,0) = x_coor(int(edge(r,1))) - x_coor(int(edge(r,0)));
        edge_vector(r,1) = y_coor(int(edge(r,1))) - y_coor(int(edge(r,0)));
    }
}

void parse_obj() {
    string obj_line;
    ifstream obj_file("/Users/garychu/Desktop/cm50245cw/ffd/ffd/man.obj");
//    ifstream obj_file("man.obj");

    while (!obj_file.eof()) {
        getline(obj_file, obj_line);
        if (obj_line[0] != 'v' & obj_line[0] != 'f') {
            continue;
        } else {
            istringstream iss(obj_line);
            vector<string> results(istream_iterator<string>{iss}, istream_iterator<string>());
            if (results[0] == "v") {
                x_coor(v_count) = stod(results[1]);
                y_coor(v_count) = -1 * stod(results[2]);
                v_count++;
            } else if (results[0] == "f") {
                f_1[f_count] = stoi(results[1]) - 1;
                f_2[f_count] = stoi(results[2]) - 1;
                f_3[f_count] = stoi(results[3]) - 1;
                f_count++;
            }
        }
    }
}

void keyboard(unsigned char key, int, int) {
    switch (key) {
        case 'q': exit(1);  break;
        //place control point
        case 'z': mouse_click_mode = 0; break;
        //select control point
        case 'x': mouse_click_mode = 1; break;
        //remove control point
        case 'c': mouse_click_mode = 2; break;
        //deform
        case 'v': mouse_click_mode = 3; break;
    }
//    glutPostRedisplay();
}

void mouse_click(int button, int state, int x, int y) {
//    double newx = (x-256.0)/256*2.5;
//    double newy = (y-256.0)/256*2.5;
    double newx = (x-256.0)/256*2;
    double newy = (y-256.0)/256*2;
    int index = find_close_point(newx, newy);
    switch (mouse_click_mode) {
        case 0: {
            //put coor into control point array/stack
            //*********use link list in c++ if can
            if (control_point[index] == -1) {
                control_point[index] = 1;
                num_control++;
            }
            break;
        }
        case 1: {
            //find close control point
            //set as handle
            if (control_point[index] == 1) {
                current_control = index;
                control_coor << x_coor(index), y_coor(index);
            }
            break;
        }
        case 2: {
            //find close control point
            //remove control point
            if (control_point[index] == 1) {
                control_point[index] = -1;
                num_control--;
                current_control = -1;
            }
            break;
        }
        case 3: {
            if (num_control > 0) {
                control_coor << newx, newy;
                sim_trans();
                calculate_T();
                scale_adj();
            }
            break;
        }
    }
    glutPostRedisplay();
}

void mouse_motion(int x, int y) {
//    double newx = (x-256.0)/256*2.5;
//    double newy = (y-256.0)/256*2.5;
//    cerr << "\t mouse is at (" << newx << ", " << newy << ")" << endl;
}

void mouse_drag(int x, int y) {
//    double newx = (x-256.0)/256*2.5;
//    double newy = (y-256.0)/256*2.5;
}

//resizing the window
void reshape(int w, int h) {
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
//    glOrtho(-2.5, 2.5, 2.5, -2.5, -1.0, 1.0);
    glOrtho(-2, 2, 2, -2, -1.0, 1.0);
    glMatrixMode(GL_MODELVIEW);
}

void display() {
    glClear(GL_COLOR_BUFFER_BIT);
    glColor3f(1.0f, 1.0f, 1.0f);

    for (int i = 0; i < 332; i++) {
        glBegin(GL_LINE_LOOP);
        glVertex2f(x_coor(f_1[i]), y_coor(f_1[i]));
        glVertex2f(x_coor(f_2[i]), y_coor(f_2[i]));
        glVertex2f(x_coor(f_3[i]), y_coor(f_3[i]));
        glEnd();
    }

    glColor3f(1.0f, 0.0f, 0.0f);
    for (int i = 0; i < 214; i++) {
        if (control_point[i] == -1) {
            continue;
        } else {
            glPointSize(5);
            glBegin(GL_POINTS);
            glVertex2f(x_coor(i), y_coor(i));
            glEnd();
        }
    }
    glutSwapBuffers();
}

void init() {
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
    parse_obj();
    find_edge();
    calculate_G();
    calculate_H();
    for (int i = 0; i < 214; i++) {
        control_point[i] = -1;
    }
}

int main(int argc, char * argv[]) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA);
    glutInitWindowSize(512, 512);
    glutInitWindowPosition(300, 100);

    glutCreateWindow("ARAP Deformation");
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutKeyboardFunc(keyboard);

    glutMouseFunc(mouse_click);
    glutPassiveMotionFunc(mouse_motion);
    glutMotionFunc(mouse_drag);

    init();
    glutMainLoop();

    return 0;
}
