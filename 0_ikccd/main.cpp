// Yifei Chen IK practice


#include <chrono>

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
    std::vector<dfm2::CVec3d> tri_pos;
    std::vector<dfm2::CVec3d> global_tripos;
    std::vector<dfm2::CVec3d> quad_pos;
    std::vector<dfm2::CVec3d> global_quadpos;
    //std::vector<dfm2::CVec3d> box_pos;
    //dfm2::CVec3d box_a;
    //dfm2::CVec3d box_b;
    //dfm2::CVec3d global_box_a;
    //dfm2::CVec3d global_box_b;
    dfm2::CVec3d globalpos;
};

dfm2::CVec3d vector_mul_matrix(dfm2::CVec3d& vec, dfm2::CMat4d& mat)
{
    dfm2::CVec3d result{ 0.0,0.0,0.0 };
    for (unsigned int i = 0; i < 3; i++)
    {
        result[i] = vec.x * mat.mat[i * 4 + 0] + vec.y * mat.mat[i * 4 + 1] + vec.z * mat.mat[i * 4 + 2];
    }
    return result;
}


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


void draw_tris(std::vector<easy_joint>& joints)
{

    ::glEnable(GL_LIGHTING);
    GLfloat cyan[] = { 0.f, .8f, .8f, 1.f };
    ::glMaterialfv(GL_FRONT, GL_DIFFUSE, cyan);
    ::glBegin(GL_TRIANGLES);
    for (unsigned int i = 0; i < joints.size(); i++)
    {
        ::glVertex3d(joints[i].global_tripos[0].x, joints[i].global_tripos[0].y, joints[i].global_tripos[0].z);
        ::glVertex3d(joints[i].global_tripos[1].x, joints[i].global_tripos[1].y, joints[i].global_tripos[1].z);
        ::glVertex3d(joints[i].global_tripos[2].x, joints[i].global_tripos[2].y, joints[i].global_tripos[2].z);
    }
    ::glEnd();
    ::glDisable(GL_LIGHTING);
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
            joints[i].globalpos = joints[i].offsetpos;
            for (unsigned int t = 0; t < 8; t++)
            {
                joints[i].global_quadpos.push_back(joints[i].globalpos + joints[i].quad_pos[t]);
            }
        }
        else
        {
            joints[i].globalpos = joints[i - 1].globalpos + joints[i].offsetpos;
            for (unsigned int t = 0; t < 8; t++)
            {
                joints[i].global_quadpos.push_back(joints[i].globalpos + joints[i].quad_pos[t]);
            }
        }
    }
}

void init_joint(std::vector<easy_joint>& joints, std::vector<double> tri_coord)
{
    std::cout << tri_coord.size() / 9 << " size" << std::endl;
    for (unsigned int i = 0; i < joints.size(); i++)
    {
        dfm2::CVec3d a{ tri_coord[i * 9 + 0],tri_coord[i * 9 + 1],tri_coord[i * 9 + 2] };
        dfm2::CVec3d b{ tri_coord[i * 9 + 3],tri_coord[i * 9 + 4],tri_coord[i * 9 + 5] };
        dfm2::CVec3d c{ tri_coord[i * 9 + 6],tri_coord[i * 9 + 7],tri_coord[i * 9 + 8] };
        joints[i].tri_pos.push_back(a);
        joints[i].tri_pos.push_back(b);
        joints[i].tri_pos.push_back(c);
        if (i == 0)
        {
            joints[i].offsetpos = joints[i].offsetpos;
            for (unsigned int t = 0; t < 3; t++)
            {
                joints[i].global_tripos.push_back(joints[i].globalpos + joints[i].tri_pos[t]);
            }
        }
        else
        {
            joints[i].offsetpos = joints[i - 1].offsetpos + joints[i].offsetpos;
            joints[i].axis = joints[i - 1].globalpos + joints[i].axis;
            for (unsigned int t = 0; t < 3; t++)
            {
                joints[i].global_tripos.push_back(joints[i].globalpos + joints[i].tri_pos[t]);
            }
        }
    }
}

void updateeasyjoint_quad(std::vector<easy_joint>& joints, int bone_index)
{
    /*
    if (bone_index == 0)
    {
        dfm2::CQuatd q0(joints[0].q);
        dfm2::CMat4d m00 = dfm2::CMat4d::Quat(q0.p);
        joints[0].globalpos = joints[0].offsetpos;
        joints[1].globalpos = joints[0].globalpos + vector_mul_matrix(joints[1].offsetpos, m00);
        for (unsigned int t = 0; t < 8; t++)
        {
            joints[0].global_quadpos[t] = vector_mul_matrix(joints[0].quad_pos[t], m00);
        }

        for (unsigned int i = 1; i < joints.size(); i++)
        {
            dfm2::CQuatd globalq(joints[0].q);
            dfm2::CQuatd global_quad_q(joints[0].q);
            dfm2::CQuatd axisq(joints[0].q);
            for (unsigned int j = 0; j <= i - 1; j++)
            {
                globalq = joints[j].q * globalq;
                if (j != i)
                {
                    axisq = joints[j].q * axisq;
                }
            }

            for (unsigned int j = 0; j <= i; j++)
            {
                global_quad_q = joints[j].q * global_quad_q;
            }
            axisq.normalize();
            globalq.normalize();
            dfm2::CMat4d m01 = dfm2::CMat4d::Quat(axisq.p);
            dfm2::CMat4d m02 = dfm2::CMat4d::Quat(globalq.p);
            dfm2::CMat4d m_global_quad = dfm2::CMat4d::Quat(global_quad_q.p);
            joints[i].globalpos = joints[i - 1].globalpos + vector_mul_matrix(joints[i].offsetpos, m02);
            joints[i].axis = joints[i - 1].globalpos + vector_mul_matrix(joints[i].axis, m01);
            joints[i].axis.normalize();

            for (unsigned int t = 0; t < 8; t++)
            {
                joints[i].global_quadpos[t] = joints[i].globalpos + vector_mul_matrix(joints[i].quad_pos[t], m_global_quad);
            }

        }
    }
    */
    //else
    //{
        for (unsigned int i = bone_index; i < joints.size(); i++)
        {
            dfm2::CQuatd globalq{1,0,0,0};
            dfm2::CQuatd global_quad_q{1,0,0,0};
            if (i != 0)
            {
                for (unsigned int j = 0; j <= i - 1; j++)
                {
                    globalq = joints[j].q * globalq;
                }
            }
                for (unsigned int z = 0; z <= i; z++)
                {
                    global_quad_q = joints[z].q * global_quad_q;
                }

                //global_quad_q = joints[i].q * global_quad_q;
                globalq.normalize();
                global_quad_q.normalize();
                //axisq.normalize();
                globalq.normalize();
                //dfm2::CMat4d m01 = dfm2::CMat4d::Quat(axisq.p);
                dfm2::CMat4d m02 = dfm2::CMat4d::Quat(globalq.p);
                dfm2::CMat4d m2_global_quad = dfm2::CMat4d::Quat(global_quad_q.p);
                if (i == bone_index)
                {
                    joints[i].globalpos = joints[i].globalpos;
                }
                else
                {
                    joints[i].globalpos = joints[i - 1].globalpos + vector_mul_matrix(joints[i].offsetpos, m02);
                }

                for (unsigned int t = 0; t < 8; t++)
                {
                    joints[i].global_quadpos[t] = joints[i].globalpos + vector_mul_matrix(joints[i].quad_pos[t], m2_global_quad);
                }
        }
}


void updateeasyjoint(std::vector<easy_joint>& joints)
{
    dfm2::CQuatd q0(joints[0].q);
    dfm2::CMat4d m01 = dfm2::CMat4d::Quat(q0.p);
    joints[0].globalpos = vector_mul_matrix(joints[0].offsetpos, m01);
    for (unsigned int t = 0; t < 3; t++)
    {
        joints[0].global_tripos[t] = joints[0].globalpos + vector_mul_matrix(joints[0].tri_pos[t], m01);
    }

    for (unsigned int i = 1; i < joints.size(); i++)
    {
        dfm2::CQuatd globalq(joints[0].q);
        for (unsigned int j = 1; j <= i; j++)
        {
            globalq = joints[j].q * globalq;
        }

        globalq.normalize();

        dfm2::CMat4d m02 = dfm2::CMat4d::Quat(globalq.p);
        joints[i].globalpos = joints[i - 1].globalpos + vector_mul_matrix(joints[i].offsetpos, m02);

        for (unsigned int t = 0; t < 3; t++)
        {
            joints[i].global_tripos[t] = joints[i - 1].globalpos + vector_mul_matrix(joints[i].tri_pos[t], m02);
        }

    }
}

void draweasybone(std::vector<easy_joint>& joints)
{


    ::glLineWidth(5);
    ::glColor3d(0.8, 0.9, 0.1);
    ::glBegin(GL_LINE_STRIP);
    for (unsigned int i = 0; i < joints.size(); i++)
    {
        ::glVertex3d(joints[i].globalpos.x, joints[i].globalpos.y, joints[i].globalpos.z);
    }
    ::glEnd();
}

void draw_bone_p(std::vector<easy_joint>& joints)
{
    ::glColor3d(0.9, 0.05, 0.05);
    for (unsigned int i = 0; i < joints.size(); i++)
    {
        ::delfem2::opengl::DrawSphereAt(32, 32, 0.5, joints[i].globalpos.x, joints[i].globalpos.y, joints[i].globalpos.z);
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

/*
void bone_init(std::vector<dfm2::CRigBone>& bones)
{
    for (unsigned int ibone = 0; ibone < bones.size(); ++ibone) {
        dfm2::CRigBone& bone = bones[ibone];
        bone.scale = 1.0;
        bone.quatRelativeRot[0] = 0.0;
        bone.quatRelativeRot[1] = 0.0;
        bone.quatRelativeRot[2] = 0.0;
        bone.quatRelativeRot[3] = 1.0;
        bone.transRelative[0] = 0.0;
        bone.transRelative[1] = 0.0;
        bone.transRelative[2] = 0.0;
        if (bone.ibone_parent != -1) {
            const dfm2::CRigBone& bone_p = bones[bone.ibone_parent];
            bone.transRelative[0] = (-bone.invBindMat[3]) - (-bone_p.invBindMat[3]);
            bone.transRelative[1] = (-bone.invBindMat[7]) - (-bone_p.invBindMat[7]);
            bone.transRelative[2] = (-bone.invBindMat[11]) - (-bone_p.invBindMat[11]);
        }
    }
    for (auto& bone : bones) {
        for (int i = 0; i < 16; ++i) { bone.affmat3Global[i] = bone.invBindMat[i]; }
        int info;
        dfm2::rig_bvh::CalcInvMat(bone.affmat3Global, 4, info);
    }
}
*/

dfm2::CQuatd qfv(dfm2::CVec3d& u, dfm2::CVec3d& v)
{
    // a concise version of q from 2 vectors
    //u.normalize();
    //v.normalize();
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
    dfm2::CVec3d target{5.0,25.0,0.0};
    
    std::vector<double> quad_local_pos{
        -1,0,1,
        -1,0,-1,
        1,0,1,
        1,0,-1,
        1,2,1,
        1,2,-1,
        -1,2,1,
        -1,2,-1,//2
        -1,0,1,
        -1,0,-1,
        1,0,1,
        1,0,-1,
        1,7,1,
        1,7,-1,
        -1,7,1,
        -1,7,-1,//3
        -1,0,1,
        -1,0,-1,
        1,0,1,
        1,0,-1,
        1,4,1,
        1,4,-1,
        -1,4,1,
        -1,4,-1,//4
        -1,0,1,
        -1,0,-1,
        1,0,1,
        1,0,-1,
        1,5,1,
        1,5,-1,
        -1,5,1,
        -1,5,-1,//5
        -1,0,1,
        -1,0,-1,
        1,0,1,
        1,0,-1,
        1,3,1,
        1,3,-1,
        -1,3,1,
        -1,3,-1,//6
        -1,0,1,
        -1,0,-1,
        1,0,1,
        1,0,-1,
        1,5,1,
        1,5,-1,
        -1,5,1,
        -1,5,-1,//7
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
    dfm2::CVec3d b1{ 0,2,0 };
    dfm2::CVec3d c1{ 0,7,0 };
    dfm2::CVec3d d1{ 0,4,0 };
    dfm2::CVec3d e1{ 0,5,0 };
    dfm2::CVec3d f1{ 0,3,0 };
    dfm2::CVec3d g1{ 0,5,0 };

    dfm2::CQuatd q{ 1,0,0,0 };
    q.normalize();

    easy_joint bone1{ "base",a1 ,q, {0,0,0},-1 };
    easy_joint bone2{ "joint1",b1 ,q,{0,0,0},0 };
    easy_joint bone3{ "joint2",c1 ,q,{0,0,0},1 };
    easy_joint bone4{ "joint3",d1 ,q,{0,0,0},2 };
    easy_joint bone5{ "joint4",e1 ,q,{0,0,0},3 };
    easy_joint bone6{ "joint5",f1 ,q,{0,0,0},4 };
    easy_joint bone7{ "end",g1 ,q,{0,0,0},5 };

    std::vector<easy_joint> all_easy_bone;
    all_easy_bone.push_back(bone1);
    all_easy_bone.push_back(bone2);
    all_easy_bone.push_back(bone3);
    all_easy_bone.push_back(bone4);
    all_easy_bone.push_back(bone5);
    all_easy_bone.push_back(bone6);
    all_easy_bone.push_back(bone7);

    dfm2::glfw::CViewer3 viewer_source(35);

    viewer_source.width = 1920;
    viewer_source.height = 1080;
    viewer_source.window_title = "IK_ccd";
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
    updateeasyjoint_quad(all_easy_bone,0);
    std::vector<double> dis_vector;
    while (!glfwWindowShouldClose(viewer_source.window))
    {

        glfwMakeContextCurrent(viewer_source.window);
        viewer_source.DrawBegin_oldGL();
        ::glEnable(GL_LINE_SMOOTH);
        ::glEnable(GL_BLEND);
        ::glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


        Floor floor{ 100, +0.1 };

        int converge_time = 1;


        for (unsigned int c = 0 ; c < converge_time; c++)
        {
            for (int i = all_easy_bone.size() - 2; i >= 0; i--)
            {
                updateeasyjoint_quad(all_easy_bone, i);
                dfm2::CVec3d jointnow(all_easy_bone[i].globalpos);
                dfm2::CVec3d end(all_easy_bone[all_easy_bone.size() - 1].globalpos);
                dfm2::CVec3d endeff_dir = end - jointnow;
                dfm2::CVec3d target_dir = target - jointnow;
                endeff_dir.normalize();
                target_dir.normalize();
                dfm2::CQuatd r_q;
                r_q = qfv(endeff_dir, target_dir);
                //r_q.SetSmallerRotation();
                r_q.normalize();

                /*
                if (i == 0)
                {
                    dfm2::CQuatd q_current = r_q * all_easy_bone[i].q;
                    q_current.normalize();
                    double y[4];
                    double xz[4];
                    SeparateYRot(y, xz, q_current.p);

                    dfm2::CQuatd Qxz(xz);
                    dfm2::CQuatd Qxz_inv = Qxz.conjugate();
                    dfm2::CQuatd q0a = q_current * Qxz_inv;
                    q0a.normalize();

                    all_easy_bone[i].q = q0a;
                }
                */
                if (i == 0)
                {
                    dfm2::CQuatd q_current = r_q * all_easy_bone[i].q;
                    q_current.normalize();
                    double y[4];
                    double xz[4];
                    if (i == 0)
                        SeparateYRot(y, xz, q_current.p);
                    else
                    {
                        SeparateXRot(y, xz, q_current.p);
                    }

                    dfm2::CQuatd Qxz(xz);
                    dfm2::CQuatd Qxz_inv = Qxz.conjugate();
                    dfm2::CQuatd q0a = q_current * Qxz_inv;
                    q0a.normalize();

                    all_easy_bone[i].q = q0a;
                }
                else if (i == 1)
                {
                    dfm2::CQuatd qnew = r_q * all_easy_bone[i].q;
                    qnew.normalize();
                    dfm2::CVec3d v_q{ qnew.x,qnew.y,qnew.z };
                    v_q = v_q / v_q.norm();
                    double theta = 2 * atan2(v_q.norm(), qnew.w);
                    double new_r = std::clamp(theta, -M_PI / 2, M_PI / 2);
                    dfm2::CQuatd qfinal{ cos(new_r / 2),v_q.x * sin(new_r / 2),v_q.y * sin(new_r / 2),v_q.z * sin(new_r / 2) };
                    qfinal.normalize();

                    all_easy_bone[i].q = qfinal;
                }
                else
                {
                    dfm2::CQuatd qnew = r_q * all_easy_bone[i].q;
                    all_easy_bone[i].q = qnew;
                }
                updateeasyjoint_quad(all_easy_bone, i);

            }
        }
        
        draw_bone_all(all_easy_bone);
        draw_box(all_easy_bone);
        draw_box_edge(all_easy_bone);

        draw_target(target);
        draw_coordinate();
        dis_vector.clear();
        for (unsigned int i = 1; i < 5; i++)
        {
            dfm2::CVec3d disv = all_easy_bone[i].globalpos - all_easy_bone[i - 1].globalpos;
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
            ImGui::Begin("IK panel");
            ImGui::Text("Control Guide");
            ImGui::Text("W S A D Q E - control the movement of target");
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