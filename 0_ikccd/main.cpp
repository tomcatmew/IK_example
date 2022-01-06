#include <chrono>

#include <iostream>
#include <vector>
#include <algorithm>
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl2.h"

#include "traj.h"
#include "draw.h"
#include "delfem2/rig_geo3.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"
#include "delfem2/opengl/old/funcs.h"
#include "delfem2/opengl/old/rigv3.h"
#include "delfem2/rig_bvh.h"

#include "math.h"
namespace dfm2 = delfem2;

struct easy_joint
{
    std::string name;
    dfm2::CVec3d pos;
    dfm2::CQuatd q;
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

void updateeasyjoint(std::vector<easy_joint>& joints)
{
    //dfm2::CMat4d m00 = dfm2::CMat4d::Quat(joints[0].q.p);
    joints[0].globalpos = joints[0].pos;
    for (unsigned int i = 1; i < joints.size(); i++)
    {
        dfm2::CQuatd globalq(joints[0].q);
        for (unsigned int j = 1; j <= i; j++)
        {
            globalq = joints[j].q * globalq;
        }


        globalq.normalize();

        dfm2::CMat4d m02 = dfm2::CMat4d::Quat(globalq.p);
        joints[i].globalpos = joints[i - 1].globalpos + vector_mul_matrix(joints[i].pos, m02);
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
    draweasybone(joints);
    draw_bone_p(joints);
}

void draw_coordinate()
{
    ::glLineWidth(5);
    ::glColor3d(0.1, 0.1, 0.98);
    ::glBegin(GL_LINE_STRIP);
    ::glVertex3d(0, 0, 0);
    ::glVertex3d(0, 0, 20);
    
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
    dfm2::CVec3d target{5.0,25.0,0.0};

    double ws_forward = 0.0f;
    double ad_leftright = 0.0f;
    double qe_topdown = 0.0f;

    //relative positions
    dfm2::CVec3d a1{0,0,0};
    dfm2::CVec3d b1{ 0,6,0 };
    dfm2::CVec3d c1{ 0,10,0 };
    dfm2::CVec3d d1{ 0,4,0 };
    dfm2::CVec3d e1{ 0,5,0 };

    dfm2::CQuatd q{ 1,0,0,0 };
    q.normalize();

    easy_joint bone1{ "base",a1 ,q };
    easy_joint bone2{ "joint1",b1 ,q };
    easy_joint bone3{ "joint2",c1 ,q };
    easy_joint bone4{ "joint3",d1 ,q };
    easy_joint bone5{ "joint4",e1 ,q };

    std::vector<easy_joint> all_easy_bone;
    all_easy_bone.push_back(bone1);
    all_easy_bone.push_back(bone2);
    all_easy_bone.push_back(bone3);
    all_easy_bone.push_back(bone4);
    all_easy_bone.push_back(bone5);

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
    updateeasyjoint(all_easy_bone);

    while (!glfwWindowShouldClose(viewer_source.window))
    {

        glfwMakeContextCurrent(viewer_source.window);
        viewer_source.DrawBegin_oldGL();
        ::glEnable(GL_LINE_SMOOTH);
        ::glEnable(GL_BLEND);
        ::glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


        Floor floor{ 100, +0.1 };

        int converge_time = 3;

        for (unsigned int j = 0; j < converge_time; j++)
        {
            for (int i = all_easy_bone.size() - 2; i > -1; i--)
            {
                dfm2::CVec3d end(all_easy_bone[all_easy_bone.size() - 1].globalpos);
                dfm2::CVec3d jointnow(all_easy_bone[i].globalpos);

                dfm2::CVec3d jointnow_dir = end - jointnow;
                dfm2::CVec3d target_dir = target - jointnow;

                dfm2::CQuatd r_q;
                r_q = qfv(jointnow_dir, target_dir);
                r_q.normalize();
                dfm2::CQuatd qnew = r_q * all_easy_bone[i + 1].q;
                //dfm2::CQuatd qnew2 = r_test2 * all_easy_bone[2].q;
                qnew.normalize();
                //qnew2.normalize();
                all_easy_bone[i + 1].q = qnew;
                updateeasyjoint(all_easy_bone);
            }
        }
        draw_bone_all(all_easy_bone);
        draw_target(target);

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

        floor.draw_checkerboard();
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
            ImGui::Separator();
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