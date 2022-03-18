// Yifei Chen IK practice jacobian transpose method


#include <chrono>
#include <filesystem>

#include <iostream>
#include <vector>
#include <algorithm>
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>
#include <chrono>
#include <thread>

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl2.h"

#include "draw.h"
#include "delfem2/rig_geo3.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/rigv3.h"
#include "delfem2/rig_bvh.h"

#include "math.h"
namespace dfm2 = delfem2;
using namespace std::this_thread; // sleep_for, sleep_until
using namespace std::chrono;
void SeparateYRot(
    double qy[4],
    double qxz[4],
    const double q[4]) {
    qxz[3] = std::sqrt(q[3] * q[3] + q[1] * q[1]);
    qy[0] = 0;
    qy[1] = q[1] / qxz[3];
    qy[2] = 0;
    qy[3] = q[3] / qxz[3];

    double det = qy[3] * qy[3] + qy[1] * qy[1];
    double invdet = 1 / det;
    qxz[0] = (qy[3] * q[0] - qy[1] * q[2]) * invdet;
    qxz[1] = 0;
    qxz[2] = (qy[1] * q[0] + qy[3] * q[2]) * invdet;


#ifndef NDEBUG
    assert(fabs(delfem2::Length_Quat(q) - 1.) < 1.0e-5);
    double qyqxz[4];
    delfem2::QuatQuat(qyqxz, qy, qxz);
    assert(fabs(q[0] - qyqxz[0]) < 1.0e-10);
    assert(fabs(q[1] - qyqxz[1]) < 1.0e-10);
    assert(fabs(q[2] - qyqxz[2]) < 1.0e-10);
    assert(fabs(q[3] - qyqxz[3]) < 1.0e-10);
#endif
}


void SeparateZRot(
    double qz[4],
    double qxy[4],
    const double q[4]) {
    qxy[3] = std::sqrt(q[3] * q[3] + q[2] * q[2]);
    qz[0] = 0;
    qz[1] = 0;
    qz[2] = q[2] / qxy[3];
    qz[3] = q[3] / qxy[3];
    //
    double det = qz[3] * qz[3] + qz[2] * qz[2];
    double invdet = 1 / det;
    qxy[0] = (qz[3] * q[0] + qz[2] * q[1]) * invdet;
    qxy[1] = (qz[3] * q[1] - qz[2] * q[0]) * invdet;
    qxy[2] = 0;
    //
#ifndef NDEBUG
    assert(fabs(dfm2::Length_Quat(q) - 1.) < 1.0e-5);
    double qyqxz[4];
    dfm2::QuatQuat(qyqxz, qz, qxy);
    assert(fabs(q[0] - qyqxz[0]) < 1.0e-10);
    assert(fabs(q[1] - qyqxz[1]) < 1.0e-10);
    assert(fabs(q[2] - qyqxz[2]) < 1.0e-10);
    assert(fabs(q[3] - qyqxz[3]) < 1.0e-10);
#endif
}

void SeparateXRot(
    double qx[4],
    double qzy[4],
    const double q[4]) {
    qzy[3] = std::sqrt(q[3] * q[3] + q[0] * q[0]);
    qx[0] = q[0] / qzy[3];
    qx[1] = 0;
    qx[2] = 0;
    qx[3] = q[3] / qzy[3];
    //
    double det = qx[3] * qx[3] + qx[0] * qx[0];
    double invdet = 1 / det;
    qzy[0] = 0;
    qzy[1] = (qx[3] * q[1] + qx[0] * q[2]) * invdet;
    qzy[2] = (qx[3] * q[2] - qx[0] * q[1]) * invdet;
    //
#ifndef NDEBUG
    assert(fabs(dfm2::Length_Quat(q) - 1.) < 1.0e-5);
    double qyqxz[4];
    dfm2::QuatQuat(qyqxz, qx, qzy);
    assert(fabs(q[0] - qyqxz[0]) < 1.0e-10);
    assert(fabs(q[1] - qyqxz[1]) < 1.0e-10);
    assert(fabs(q[2] - qyqxz[2]) < 1.0e-10);
    assert(fabs(q[3] - qyqxz[3]) < 1.0e-10);
#endif
}

void draw_arrow(dfm2::CVec3d& p0, dfm2::CVec3d& p1)
{
    ::glEnable(GL_LIGHTING);
    delfem2::opengl::DrawArrowOcta_FaceNrm(
        p0,
        p1,
        0.1,
        0.2);
    //
    ::glDisable(GL_LIGHTING);
}

void draw_cube(double a[], double b[]) {
    //::glColor3d(0.34, 0.45, 0.61);
    ::glEnable(GL_LIGHTING);
    GLfloat cyan[] = { 0.f, .8f, .8f, 1.f };
    ::glMaterialfv(GL_FRONT, GL_DIFFUSE, cyan);
    dfm2::opengl::DrawBox3_Face(a, b);
    ::glDisable(GL_LIGHTING);
}

struct easy_joint
{
    std::string name;
    dfm2::CVec3d offsetpos;
    dfm2::CQuatd q;
    dfm2::CVec3d axis;
    int iparent;
    double affineglobal[16] ={1,0,0,0,
                              0,1,0,0,
                              0,0,1,0,
                              0,0,0,1};
    std::vector<dfm2::CVec3d> quad_pos;
    std::vector<dfm2::CVec3d> global_quadpos;

    dfm2::CVec3d globalpos()
    {
        return dfm2::CVec3d(affineglobal[3], affineglobal[7], affineglobal[11]);
    }
        
};


void draw_box(std::vector<easy_joint>& joints) {
    ::glEnable(GL_LIGHTING);
    for (unsigned int i = 0; i < joints.size(); i++)
    {
        GLfloat cyan[] = { 0.f, .8f, .8f, 1.f };
        GLfloat cred[] = { 0.9f, .1f, .12f, 1.f };
        GLfloat cgreen[] = { 0.5f, 0.5f, .1f, 1.f };

        if (i == 0)
        {
            ::glMaterialfv(GL_FRONT, GL_DIFFUSE, cgreen);
        }
        else if (i == 1)
        {
            ::glMaterialfv(GL_FRONT, GL_DIFFUSE, cred);
        }
        else
        {
            ::glMaterialfv(GL_FRONT, GL_DIFFUSE, cyan);
        }
        ::glBegin(GL_QUADS);
        ::glVertex3d(joints[i].global_quadpos[0].x, joints[i].global_quadpos[0].y, joints[i].global_quadpos[0].z);
        ::glVertex3d(joints[i].global_quadpos[1].x, joints[i].global_quadpos[1].y, joints[i].global_quadpos[1].z);
        ::glVertex3d(joints[i].global_quadpos[3].x, joints[i].global_quadpos[3].y, joints[i].global_quadpos[3].z);
        ::glVertex3d(joints[i].global_quadpos[2].x, joints[i].global_quadpos[2].y, joints[i].global_quadpos[2].z); //1

        ::glVertex3d(joints[i].global_quadpos[2].x, joints[i].global_quadpos[2].y, joints[i].global_quadpos[2].z);
        ::glVertex3d(joints[i].global_quadpos[4].x, joints[i].global_quadpos[4].y, joints[i].global_quadpos[4].z);
        ::glVertex3d(joints[i].global_quadpos[5].x, joints[i].global_quadpos[5].y, joints[i].global_quadpos[5].z);
        ::glVertex3d(joints[i].global_quadpos[3].x, joints[i].global_quadpos[3].y, joints[i].global_quadpos[3].z); //2

        ::glVertex3d(joints[i].global_quadpos[4].x, joints[i].global_quadpos[4].y, joints[i].global_quadpos[4].z);
        ::glVertex3d(joints[i].global_quadpos[6].x, joints[i].global_quadpos[6].y, joints[i].global_quadpos[6].z);
        ::glVertex3d(joints[i].global_quadpos[7].x, joints[i].global_quadpos[7].y, joints[i].global_quadpos[7].z);
        ::glVertex3d(joints[i].global_quadpos[5].x, joints[i].global_quadpos[5].y, joints[i].global_quadpos[5].z); //3

        ::glVertex3d(joints[i].global_quadpos[6].x, joints[i].global_quadpos[6].y, joints[i].global_quadpos[6].z);
        ::glVertex3d(joints[i].global_quadpos[0].x, joints[i].global_quadpos[0].y, joints[i].global_quadpos[0].z);
        ::glVertex3d(joints[i].global_quadpos[1].x, joints[i].global_quadpos[1].y, joints[i].global_quadpos[1].z);
        ::glVertex3d(joints[i].global_quadpos[7].x, joints[i].global_quadpos[7].y, joints[i].global_quadpos[7].z); //4

        ::glVertex3d(joints[i].global_quadpos[0].x, joints[i].global_quadpos[0].y, joints[i].global_quadpos[0].z);
        ::glVertex3d(joints[i].global_quadpos[6].x, joints[i].global_quadpos[6].y, joints[i].global_quadpos[6].z);
        ::glVertex3d(joints[i].global_quadpos[4].x, joints[i].global_quadpos[4].y, joints[i].global_quadpos[4].z);
        ::glVertex3d(joints[i].global_quadpos[2].x, joints[i].global_quadpos[2].y, joints[i].global_quadpos[2].z); //5 

        ::glVertex3d(joints[i].global_quadpos[3].x, joints[i].global_quadpos[3].y, joints[i].global_quadpos[3].z);
        ::glVertex3d(joints[i].global_quadpos[1].x, joints[i].global_quadpos[1].y, joints[i].global_quadpos[1].z);
        ::glVertex3d(joints[i].global_quadpos[7].x, joints[i].global_quadpos[7].y, joints[i].global_quadpos[7].z);
        ::glVertex3d(joints[i].global_quadpos[5].x, joints[i].global_quadpos[5].y, joints[i].global_quadpos[5].z); //6
        ::glEnd();
    }
    
    ::glDisable(GL_LIGHTING);
}

void draw_box_edge(std::vector<easy_joint>& joints) {  
    for (unsigned int i = 0; i < joints.size(); i++)
    {
        ::glLineWidth(1);
        ::glColor3d(0.1, 0.1, 0.1);
        ::glBegin(GL_LINE_STRIP);
        ::glVertex3d(joints[i].global_quadpos[0].x, joints[i].global_quadpos[0].y, joints[i].global_quadpos[0].z);
        ::glVertex3d(joints[i].global_quadpos[1].x, joints[i].global_quadpos[1].y, joints[i].global_quadpos[1].z);
        ::glVertex3d(joints[i].global_quadpos[3].x, joints[i].global_quadpos[3].y, joints[i].global_quadpos[3].z);
        ::glVertex3d(joints[i].global_quadpos[2].x, joints[i].global_quadpos[2].y, joints[i].global_quadpos[2].z);

        ::glVertex3d(joints[i].global_quadpos[0].x, joints[i].global_quadpos[0].y, joints[i].global_quadpos[0].z);
        ::glVertex3d(joints[i].global_quadpos[6].x, joints[i].global_quadpos[6].y, joints[i].global_quadpos[6].z);
        ::glVertex3d(joints[i].global_quadpos[4].x, joints[i].global_quadpos[4].y, joints[i].global_quadpos[4].z);
        ::glVertex3d(joints[i].global_quadpos[2].x, joints[i].global_quadpos[2].y, joints[i].global_quadpos[2].z);

        ::glVertex3d(joints[i].global_quadpos[3].x, joints[i].global_quadpos[3].y, joints[i].global_quadpos[3].z);
        ::glVertex3d(joints[i].global_quadpos[5].x, joints[i].global_quadpos[5].y, joints[i].global_quadpos[5].z);
        ::glVertex3d(joints[i].global_quadpos[4].x, joints[i].global_quadpos[4].y, joints[i].global_quadpos[4].z);
        ::glVertex3d(joints[i].global_quadpos[6].x, joints[i].global_quadpos[6].y, joints[i].global_quadpos[6].z);

        ::glVertex3d(joints[i].global_quadpos[7].x, joints[i].global_quadpos[7].y, joints[i].global_quadpos[7].z);
        ::glVertex3d(joints[i].global_quadpos[5].x, joints[i].global_quadpos[5].y, joints[i].global_quadpos[5].z);
        ::glVertex3d(joints[i].global_quadpos[3].x, joints[i].global_quadpos[3].y, joints[i].global_quadpos[3].z);
        ::glVertex3d(joints[i].global_quadpos[1].x, joints[i].global_quadpos[1].y, joints[i].global_quadpos[1].z);

        ::glVertex3d(joints[i].global_quadpos[7].x, joints[i].global_quadpos[7].y, joints[i].global_quadpos[7].z);
        ::glVertex3d(joints[i].global_quadpos[6].x, joints[i].global_quadpos[6].y, joints[i].global_quadpos[6].z);
        ::glVertex3d(joints[i].global_quadpos[0].x, joints[i].global_quadpos[0].y, joints[i].global_quadpos[0].z);
        ::glEnd();
    }
  
}

void init_joint_quad(std::vector<easy_joint>& joints, std::vector<double> quad_coord)
{
    std::cout << quad_coord.size() / 24 << " - size" << std::endl;
    for (unsigned int i = 0; i < joints.size(); i++)
    {
        dfm2::CVec3d a{ quad_coord[i * 24 + 0],quad_coord[i * 24 + 1],quad_coord[i * 24 + 2] };
        dfm2::CVec3d b{ quad_coord[i * 24 + 3],quad_coord[i * 24 + 4],quad_coord[i * 24 + 5] };
        dfm2::CVec3d c{ quad_coord[i * 24 + 6],quad_coord[i * 24 + 7],quad_coord[i * 24 + 8] };
        dfm2::CVec3d d{ quad_coord[i * 24 + 9],quad_coord[i * 24 + 10],quad_coord[i * 24 + 11] };
        dfm2::CVec3d e{ quad_coord[i * 24 + 12],quad_coord[i * 24 + 13],quad_coord[i * 24 + 14] };
        dfm2::CVec3d f{ quad_coord[i * 24 + 15],quad_coord[i * 24 + 16],quad_coord[i * 24 + 17] };
        dfm2::CVec3d g{ quad_coord[i * 24 + 18],quad_coord[i * 24 + 19],quad_coord[i * 24 + 20] };
        dfm2::CVec3d h{ quad_coord[i * 24 + 21],quad_coord[i * 24 + 22],quad_coord[i * 24 + 23] };
        joints[i].quad_pos.push_back(a);
        joints[i].quad_pos.push_back(b);
        joints[i].quad_pos.push_back(c);
        joints[i].quad_pos.push_back(d);
        joints[i].quad_pos.push_back(e);
        joints[i].quad_pos.push_back(f);
        joints[i].quad_pos.push_back(g);
        joints[i].quad_pos.push_back(h);
   
        if (i == 0)
        {
            //joints[i].globalpos = joints[i].offsetpos;
            for (unsigned int t = 0; t < 8; t++)
            {
                joints[i].global_quadpos.push_back(joints[i].globalpos() + joints[i].quad_pos[t]);
            }
        }
        else
        {
            //joints[i].globalpos = joints[i - 1].globalpos + joints[i].offsetpos;
            for (unsigned int t = 0; t < 8; t++)
            {
                joints[i].global_quadpos.push_back(joints[i].globalpos() + joints[i].quad_pos[t]);
            }
        }
    }
}

void update_j_final(std::vector<easy_joint>& joints, int bone_index)
{
    for (unsigned int i = bone_index; i < joints.size(); i ++)
    {
        dfm2::CMat4d m01 = dfm2::CMat4d::Translation({ joints[i].offsetpos.x,joints[i].offsetpos.y,joints[i].offsetpos.z });
        m01 = m01 * dfm2::CMat4d::Quat(joints[i].q.p);

        const int ibone_p = joints[i].iparent;
        if (ibone_p < 0 || ibone_p >= (int)joints.size())
        {
            dfm2::Copy_Mat4(joints[i].affineglobal, m01.mat);
            continue;
        }
        //assert(ibone_p < (int)ibone);

        dfm2::MatMat4(joints[i].affineglobal, joints[ibone_p].affineglobal, m01.mat);
    }
    for (unsigned int i = bone_index; i < joints.size(); i++)
    {
        for (unsigned int t = 0; t < 8; t++)
        {
            double quad_off[4] = { joints[i].quad_pos[t].x,joints[i].quad_pos[t].y,joints[i].quad_pos[t].z,1 };
            double quad_new[4] = { 0,0,0,0 };
            dfm2::MatVec4(quad_new, joints[i].affineglobal, quad_off);
            //dfm2::CVec3d quad_glo_new(quad_new[0], quad_new[1], quad_new[2]);
            joints[i].global_quadpos[t] = { quad_new[0], quad_new[1], quad_new[2] };
        }
    }
}

void reset_pose(std::vector<easy_joint>& joints)
{
    dfm2::CQuatd oq{ 1,0,0,0 };
    for (unsigned int i = 0; i < joints.size(); i++)
    {
        joints[i].q = oq;
    }
    update_j_final(joints, 0);
}

void draweasybone(std::vector<easy_joint>& joints)
{


    ::glLineWidth(5);
    ::glColor3d(0.8, 0.9, 0.1);
    ::glBegin(GL_LINE_STRIP);
    for (unsigned int i = 0; i < joints.size(); i++)
    {
        ::glVertex3d(joints[i].globalpos().x, joints[i].globalpos().y, joints[i].globalpos().z);
    }
    ::glEnd();
}

void draw_bone_p(std::vector<easy_joint>& joints)
{
    ::glColor3d(0.9, 0.05, 0.05);
    for (unsigned int i = 0; i < joints.size(); i++)
    {
        ::delfem2::opengl::DrawSphereAt(32, 32, 0.5, joints[i].globalpos().x, joints[i].globalpos().y, joints[i].globalpos().z);
    }
    
}

void draw_target(dfm2::CVec3d& po)
{
    ::glColor3d(0.1, 0.05, 0.08);
    ::delfem2::opengl::DrawSphereAt(32, 32, 0.5, po.x, po.y, po.z);
}

void draw_bone_all(std::vector<easy_joint>& joints)
{
    //draweasybone(joints);
    draw_bone_p(joints);
}

void draw_coordinate()
{
    ::glLineWidth(3);
    ::glColor3d(0.1, 0.1, 0.98);
    ::glBegin(GL_LINE_STRIP);
    ::glVertex3d(0, 0, 0);
    ::glVertex3d(0, 0, 20);
    ::glEnd();

    ::glColor3d(0.9, 0.1, 0.11);
    ::glBegin(GL_LINE_STRIP);
    ::glVertex3d(0, 0, 0);
    ::glVertex3d(20, 0, 0);
    ::glEnd();

    ::glColor3d(0.4, 0.4, 0.2);
    ::glBegin(GL_LINE_STRIP);
    ::glVertex3d(0, 0, 0);
    ::glVertex3d(0, 20, 0);
    ::glEnd();
}



void joystick_callback(int jid, int event)
{
    if (event == GLFW_CONNECTED)
    {
        std::cout << "Controller Connected" << std::endl;
        // The joystick was connected
    }
    else if (event == GLFW_DISCONNECTED)
    {
        std::cout << "Controller Disconnected" << std::endl;
        std::cout << "Recommend connecting a controller for play" << std::endl;

        // The joystick was disconnected
    }
}

static void error_callback(int error, const char* description) {
    fputs(description, stderr);
}

static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
        glfwSetWindowShouldClose(window, GL_TRUE);
}

dfm2::CQuatd qfv(dfm2::CVec3d& u, dfm2::CVec3d& v)
{
    // a concise version of q from 2 vectors
    u.normalize();
    v.normalize();
    double k = 1.0 + u.dot(v);
    double s = 1 / sqrt(k + k);
    dfm2::CVec3d cross = u.cross(v);
    cross = s * cross;
    return dfm2::CQuatd(k * s, cross.x, cross.y, cross.z);
}


dfm2::CQuatd quatfromtwovectors(dfm2::CVec3d& u, dfm2::CVec3d& v)
{
    u.normalize();
    v.normalize();
    float norm_u_norm_v = sqrt(u.norm() * v.norm());
    float real_part = norm_u_norm_v + u.dot(v);

    dfm2::CVec3d w;
    if (real_part < 1.e-6f * norm_u_norm_v) {
        /* If u and v are exactly opposite, rotate 180 degrees
         * around an arbitrary orthogonal axis. Axis normalisation
         * can happen later, when we normalise the quaternion. */
        real_part = 0.0f;
        w = abs(u.x) > abs(u.z) ? dfm2::CVec3d(-u.y, u.x, 0.f)
            : dfm2::CVec3d(0.f, -u.z, u.y);
    }
    else {
        /* Otherwise, build quaternion the standard way. */
        w = u.cross(v);
    }

    return dfm2::CQuatd(real_part, w.x, w.y, w.z);
}


dfm2::CQuatd old_quatfromtwovectors(dfm2::CVec3d& u, dfm2::CVec3d& v)
{
    u.normalize();
    v.normalize();
    double dot = u.dot(v);
    dfm2::CVec3d axis = u.cross(v);
    double axis_norm = axis.norm();
    axis = axis / axis_norm;
    //printf("x: %f\n", dot);
    //printf("y: %f\n", axis_norm);
    double angle_a = atan(axis_norm / dot);
    //double angle_a = atan2(axis_norm, dot);
    //printf("Angle: %f\n", angle_a);
    dfm2::CQuatd result{ cos(angle_a / 2), axis.x * sin(angle_a / 2) ,axis.y * sin(angle_a / 2) ,axis.z * sin(angle_a / 2) };
    return result;

}





// =======================================
int main(int argc, char* argv[]) {
    dfm2::CVec3d target{5.0,30.0,0.0};
    
    std::vector<double> quad_local_pos{
        -1,0,1,
        -1,0,-1,
        1,0,1,
        1,0,-1,
        1,8,1,
        1,8,-1,
        -1,8,1,
        -1,8,-1,//1
        -1,0,1,
        -1,0,-1,
        1,0,1,
        1,0,-1,
        1,8,1,
        1,8,-1,
        -1,8,1,
        -1,8,-1,//2
        -1,0,1,
        -1,0,-1,
        1,0,1,
        1,0,-1,
        1,8,1,
        1,8,-1,
        -1,8,1,
        -1,8,-1, //3
        0,0,0,
        0,0,0,
        0,0,0,
        0,0,0,
        0,0,0,
        0,0,0,
        0,0,0,
        0,0,0  //endeffector
    };
   
    double ws_forward = 0.0f;
    double ad_leftright = 0.0f;
    double qe_topdown = 0.0f;

    //relative positions
    dfm2::CVec3d a1{0,0,0};
    dfm2::CVec3d b1{ 0,8,0 };
    dfm2::CVec3d c1{ 0,8,0 };
    dfm2::CVec3d d1{ 0,8,0 };

    dfm2::CQuatd q{ 1,0,0,0 };
    q.normalize();

    easy_joint bone1{ "base",a1 ,q, {0,0,0},-1 };
    easy_joint bone2{ "joint1",b1 ,q,{0,0,0},0 };
    easy_joint bone3{ "joint2",c1 ,q,{0,0,0},1 };
    easy_joint bone4{ "end",d1 ,q,{0,0,0},2 };

    std::vector<easy_joint> all_easy_bone;
    all_easy_bone.push_back(bone1);
    all_easy_bone.push_back(bone2);
    all_easy_bone.push_back(bone3);
    all_easy_bone.push_back(bone4);

    dfm2::glfw::CViewer3 viewer_source(35);

    viewer_source.width = 1920;
    viewer_source.height = 1080;
    viewer_source.window_title = "IK_Jacobian";
    viewer_source.view_rotation = std::make_unique<delfem2::ModelView_Ytop>();

    delfem2::glfw::InitGLOld();

    viewer_source.OpenWindow();

    delfem2::opengl::setSomeLighting();

    std::cout << "run starts" << std::endl;

    glfwSetErrorCallback(error_callback);
    if (!glfwInit())
        exit(EXIT_FAILURE);

    if (!viewer_source.window)
    {
        glfwTerminate();
        exit(EXIT_FAILURE);
    }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);

    glfwMakeContextCurrent(viewer_source.window);

    //======
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    io.Fonts->AddFontDefault();
    //std::filesystem::path path = std::filesystem::current_path();
    //std::string fontpath = path.string() + "/../Roboto-Medium.ttf";
    //io.Fonts->AddFontFromFileTTF(fontpath.c_str(), 20.0f);
    //io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;     // Enable Keyboard Controls
    //io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;      // Enable Gamepad Controls

    // Setup Dear ImGui style
    ImGui::StyleColorsDark();
    //ImGui::StyleColorsClassic();

    // Setup Platform/Renderer backends
    ImGui_ImplGlfw_InitForOpenGL(viewer_source.window, true);
    ImGui_ImplOpenGL2_Init();

    // Our state
    bool show_demo_window = true;
    bool show_another_window = false;
    ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);
    //=======
    //assert((tri_local_pos.size() == 451) && "Yet another way to add assert message");
    init_joint_quad(all_easy_bone, quad_local_pos);
    //init_joint(all_easy_bone, tri_local_pos);
    //updateeasyjoint_quad(all_easy_bone,0);
    update_j_final(all_easy_bone,0);
    std::vector<double> dis_vector;

    dfm2::CMat3d jaco = dfm2::CMat3d::Identity();
    jaco.p_[0] = 9.9f;
    for (unsigned int i = 0; i < 9; i++)
    {
        std::cout << jaco.p_[i] << std::endl;
    }

    while (!glfwWindowShouldClose(viewer_source.window))
    {

        glfwMakeContextCurrent(viewer_source.window);
        viewer_source.DrawBegin_oldGL();
        ::glEnable(GL_LINE_SMOOTH);
        ::glEnable(GL_BLEND);
        ::glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


        Floor floor{ 100, +0.1 };

        int converge_time = 30;

        //reset_pose(all_easy_bone);

        // Jacobian Transpose
        //dfm2::CVec3d end(all_easy_bone[all_easy_bone.size() - 1].globalpos());
        //double h = 0.00001;


        double beta = 0.001;
        dfm2::CVec3d old_end(all_easy_bone[all_easy_bone.size() - 1].globalpos());
        dfm2::CVec3d distance_eg = target - old_end;

        unsigned int count = 0;
        //  currently I let it iterate 200 times, I should check the distance and stop the loop when target is close enough
        // but I failed to use the distance stop the loop, I may mess up some parts. 
        while (count < 50)
        {
            // iterate through 3 axis, Z, Y ,X
            // currently only works for Z axis, for 3 axis I need to change axis based on the parent's bone, I haven't finish yet 
            for (unsigned int i = 0; i < 1; i++)
            {
                //get end effector pos
                dfm2::CVec3d end(all_easy_bone[all_easy_bone.size() - 1].globalpos());
                // compute delta e 
                distance_eg = target - end;
                // mutiply a small factor beta
                dfm2::CVec3d dis_end_target = beta * distance_eg;
                dfm2::CVec3d axis;
              
                // iterate through 3 axis, Z, Y ,X
                if (i == 0)
                    axis = { 0.0,0.0,1.0 };
                else if (i == 1)
                    axis = { 0.0,1.0,0.0 };
                else
                    axis = { 1.0,0.0,0.0 };

                //position of bone 0, 1, 2
                dfm2::CVec3d p0 = all_easy_bone[0].globalpos();
                dfm2::CVec3d p1 = all_easy_bone[1].globalpos();
                dfm2::CVec3d p2 = all_easy_bone[2].globalpos();

                //vector between end effector to bone 0,1,2
                dfm2::CVec3d v0 = end - p0;
                dfm2::CVec3d v1 = end - p1;
                dfm2::CVec3d v2 = end - p2;

                // create jacobian matrix
                dfm2::CMat3d jaco = dfm2::CMat3d::Identity();

                // using analytical method to get jacobian matrix
                dfm2::CVec3d J_A = axis.cross(v0);
                dfm2::CVec3d J_B = axis.cross(v1);
                dfm2::CVec3d J_C = axis.cross(v2);

                jaco.p_[0] = J_A.x;
                jaco.p_[1] = J_B.x;
                jaco.p_[2] = J_C.x;

                jaco.p_[3] = J_A.y;
                jaco.p_[4] = J_B.y;
                jaco.p_[5] = J_C.y;

                jaco.p_[6] = J_A.z;
                jaco.p_[7] = J_B.z;
                jaco.p_[8] = J_C.z;

                // get jacobian transpose
                dfm2::CMat3d jaco_t = jaco.transpose();
                dfm2::CVec3d d_theta;
                // mutiple jacobian tranpose with delta e, compute change in joint DOF
                dfm2::MatVec3(d_theta.p, jaco_t.p_, dis_end_target.p);

                d_theta = d_theta * 0.05;
                // Geting the change in rotation, delta_rv0,rv1,rv2
                // mutiply it to the rotation axis, so we have rotation vector with magnitude (vector3 of axis-angle form)
                dfm2::CVec3d delta_rv0 = axis * d_theta.x;
                dfm2::CVec3d delta_rv1 = axis * d_theta.y;
                dfm2::CVec3d delta_rv2 = axis * d_theta.z;

                // convert the axis-angle form to quaternion
                dfm2::CQuatd delta_r0;
                //dfm2::Quaternion_EulerAngle(delta_r0.p, {0.f,0.f,d_theta.x}, { 0,1,2 });
                dfm2::Quat_CartesianAngle(delta_r0.p, delta_rv0.p);
                dfm2::CQuatd delta_r1;
                //dfm2::Quaternion_EulerAngle(delta_r1.p, { 0.f,0.f,d_theta.y }, { 0,1,2 });
                dfm2::Quat_CartesianAngle(delta_r1.p, delta_rv1.p);
                dfm2::CQuatd delta_r2;
                //dfm2::Quaternion_EulerAngle(delta_r2.p, { 0.f,0.f,d_theta.z }, { 0,1,2 });
                dfm2::Quat_CartesianAngle(delta_r2.p, delta_rv2.p);
                
                // mutiply quaternion to bone0,1,2 local rotation quaternion, then update their global positions
                
                dfm2::CQuatd qnew0 = delta_r0 * all_easy_bone[0].q;
                all_easy_bone[0].q = qnew0;
                update_j_final(all_easy_bone, 0);

                
                dfm2::CQuatd qnew1 = delta_r1 * all_easy_bone[1].q;
                all_easy_bone[1].q = qnew1;
                update_j_final(all_easy_bone, 1);

                
                dfm2::CQuatd qnew3 = delta_r2 * all_easy_bone[2].q;
                all_easy_bone[2].q = qnew3;
                update_j_final(all_easy_bone, 2);

                //update end effector
                update_j_final(all_easy_bone, 3);
            }
            count++;
        }

        draw_bone_all(all_easy_bone);
        draw_box(all_easy_bone);
        draw_box_edge(all_easy_bone);

        draw_target(target);
        draw_coordinate();
        dis_vector.clear();
        for (unsigned int i = 1; i < 5; i++)
        {
            dfm2::CVec3d disv = all_easy_bone[i].globalpos() - all_easy_bone[i - 1].globalpos();
            double dis = disv.norm();
            dis_vector.push_back(dis);
        }

        double a[3] = { 0,0,0 };
        double b[3] = { 1,1,1 };
        //draw_cube(a, b);
        // 
        //dfm2::opengl::DrawBone_Octahedron
        // camera and gamepad controls

        float zoom_speed = 0.02;
        int present = 0;
        float drop_speed = 0.2;
        float gamepad_x = 0.0f;
        float gamepad_y = 0.0f;
        float gamepad2_x = 0.0f;
        float gamepad2_y = 0.0f;
        bool buttonLB = false;
        bool buttonRB = false;
        bool check_controller = false;

       

        glfwSetJoystickCallback(joystick_callback);
        if (glfwJoystickPresent(GLFW_JOYSTICK_1) == GLFW_TRUE)
        {

            present = glfwJoystickPresent(GLFW_JOYSTICK_1);

            if (present == 1)
            {
                check_controller = true;
                int axesCount;
                const float* axes = glfwGetJoystickAxes(GLFW_JOYSTICK_1, &axesCount);
                gamepad_x = axes[0];
                gamepad_y = axes[1];
                gamepad2_x = abs(axes[2]) < 0.05 ? 0.0f : axes[2] / 100;
                gamepad2_y = abs(axes[3]) < 0.05 ? 0.0f : axes[3] / 100;
                float gamepad_mag = sqrtf(gamepad_x * gamepad_x + gamepad_y * gamepad_y);

                int triggerRT = axes[5];
                int triggerLT = axes[4];
                if (gamepad_mag > 0.2f)
                {
                    float gamepaddirx = gamepad_x / gamepad_mag;
                    float gamepaddiry = gamepad_y / gamepad_mag;

                    target.x += gamepaddirx;
                    target.z += gamepaddiry;

                }
                if (triggerLT == 1) {
                    target.y += -drop_speed;
                }
                if (triggerRT == 1) {
                    target.y += drop_speed;
                }

                int buttonCont;
                const unsigned char* buttons = glfwGetJoystickButtons(GLFW_JOYSTICK_1, &buttonCont);
                //0 = A, 1 = B, 2 = X, 3 = Y
                if (GLFW_PRESS == buttons[4])
                {
                    buttonLB = true;
                }
                else if (GLFW_RELEASE == buttons[4])
                {
                    buttonLB = false;
                }
                if (GLFW_PRESS == buttons[5])
                {
                    buttonRB = true;
                }
                else if (GLFW_RELEASE == buttons[5])
                {
                    buttonRB = false;
                }

                if (buttonRB == true)
                {
                    viewer_source.scale += zoom_speed;
                }
                if (buttonLB == true)
                {
                    viewer_source.scale -= zoom_speed;
                }
            }

        }
        else
        {
            int state_w = glfwGetKey(viewer_source.window, GLFW_KEY_W);
            int state_s = glfwGetKey(viewer_source.window, GLFW_KEY_S);
            int state_a = glfwGetKey(viewer_source.window, GLFW_KEY_A);
            int state_d = glfwGetKey(viewer_source.window, GLFW_KEY_D);
            int state_q = glfwGetKey(viewer_source.window, GLFW_KEY_Q);
            int state_e = glfwGetKey(viewer_source.window, GLFW_KEY_E);

            //now only support moving on YX plane
            /*
            if (state_w == GLFW_PRESS)
            {
                ws_forward = -0.3f;
            }
            if (state_s == GLFW_PRESS)
            {
                ws_forward = 0.3f;
            }
            if ((state_s == GLFW_RELEASE) && (state_w == GLFW_RELEASE))
            {
                ws_forward = 0.0f;
            }
            */
            if (state_a == GLFW_PRESS)
            {
                ad_leftright = -0.3f;
            }
            if (state_d == GLFW_PRESS)
            {
                ad_leftright = 0.3f;
            }
            if ((state_d == GLFW_RELEASE) && (state_a == GLFW_RELEASE))
            {
                ad_leftright = 0.0f;
            }
            if (state_q == GLFW_PRESS)
            {
                qe_topdown = 0.3f;
            }
            if (state_e == GLFW_PRESS)
            {
                qe_topdown = -0.3f;
            }
            if ((state_e == GLFW_RELEASE) && (state_q == GLFW_RELEASE))
            {
                qe_topdown = 0.0f;
            }
        }

        target.x += ad_leftright;
        target.y += qe_topdown;
        target.z += ws_forward;

       
        viewer_source.view_rotation->Rot_Camera(-static_cast<float>(gamepad2_x), -static_cast<float>(gamepad2_y));

        //floor.draw_checkerboard();
        // GUI *******
        ImGui_ImplOpenGL2_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();
        // 1. Show the big demo window (Most of the sample code is in ImGui::ShowDemoWindow()! You can browse its code to learn more about Dear ImGui!).
        {
            //static float f = 0.0f;
            ImGui::Begin("IK Jacobian panel");
            ImGui::Text("Control Guide");
            ImGui::Text("A D Q E - control the movement of target");
            ImGui::Text("You cannot use W S because it only allows rotation on Z aixs.");
            ImGui::Text("I'm currently working on to make it support X and Y axis rotation");
            ImGui::NewLine();
            ImGui::Separator();
            ImGui::NewLine();
            ImGui::Text("Length of Bone - 0: %f ; 1: %f ; 3: %f ; 4: %f", dis_vector[0], dis_vector[1], dis_vector[2], dis_vector[3]);
            ImGui::NewLine();
            ImGui::Text("Target Position - X: %f ; Y : %f, Z: %f", target.x, target.y, target.z);
            ImGui::Separator();
            ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
            ImGui::End();
        }

        // Rendering
        ImGui::Render();
        int display_w, display_h;
        glfwGetFramebufferSize(viewer_source.window, &display_w, &display_h);
        glViewport(0, 0, display_w, display_h);
        //glClearColor(clear_color.x * clear_color.w, clear_color.y * clear_color.w, clear_color.z * clear_color.w, clear_color.w);
        //glClear(GL_COLOR_BUFFER_BIT);
        ImGui_ImplOpenGL2_RenderDrawData(ImGui::GetDrawData());
        // GUI  *********
        viewer_source.SwapBuffers();
        //glfwSwapBuffers(viewer_source.window);
        glfwPollEvents();
    }


    glfwDestroyWindow(viewer_source.window);
    glfwTerminate();
    exit(EXIT_SUCCESS);
}