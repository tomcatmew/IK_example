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

namespace dfm2 = delfem2;

void draw_coordinate()
{
    ::glLineWidth(5);
    ::glColor3d(0.1, 0.1, 0.98);
    ::glBegin(GL_LINE_STRIP);
    ::glVertex3d(0, 0, 0);
    ::glVertex3d(0, 0, 20);
    
    ::glEnd();
}

void unrotate_Y_future_traj(std::vector<double>& future_traj, dfm2::CVec2d& face_dir, std::vector<double>& unY_future_traj)
{

    dfm2::CQuatd q0y;
    const dfm2::CMat3d R0(
        face_dir.y, 0, face_dir.x,
        0, 1, 0,
        -face_dir.x, 0, face_dir.y);
    q0y = R0.GetQuaternion();
    q0y.SetSmallerRotation();
    dfm2::CQuatd q0yinv = q0y.conjugate();
    for (unsigned int i = 0; i < future_traj.size() / 2; i++)
    {
        double pos[3] = { future_traj[i * 2 + 0],0.0,future_traj[i * 2 + 1] };
        double new_pos[3];
        dfm2::QuatVec(new_pos, q0yinv.p, pos);
        unY_future_traj.push_back(new_pos[0]);
        unY_future_traj.push_back(new_pos[2]);
    }
}

void debug_future_traj_correct(std::vector<double>& traj_database, dfm2::CVec2d& face_dir)
{
    ::glLineWidth(4);
    ::glColor3d(1.0, 1.0, 0);
    ::glBegin(GL_LINE_STRIP);
        ::glVertex3d(0.0, 0.1, 0.0);
        ::glVertex3d(face_dir.x, 0.1, face_dir.y);
    ::glEnd();

    for (int j = 0; j < 5; j++)
    {
        float color[3] = { 0.98f,0.1f,0.1f };
        ::glColor3fv(color);
        ::delfem2::opengl::DrawSphereAt(64, 64, 0.7, traj_database[j * 2 + 0], 0, traj_database[j * 2 + 1]);
    }

}

void debug_traj_correct(std::vector<double>& traj_database, int nframe)
{

    for (int j = 0; j < 5; j++)
    {
        float color[3] = { 0.1f,0.1f,0.1f };
        ::glColor3fv(color);
        ::delfem2::opengl::DrawSphereAt(64, 64, 0.7, traj_database[nframe * 10 + j * 2 + 0], 0, traj_database[nframe * 10 + j * 2 + 1]);
    }
    
}

void SeparateYRot(
    double qy[4],
    double qxz[4],
    const double q[4]) {
    qxz[3] = std::sqrt(q[3] * q[3] + q[1] * q[1]);
    qy[0] = 0;
    qy[1] = q[1] / qxz[3];
    qy[2] = 0;
    qy[3] = q[3] / qxz[3];
    //
    double det = qy[3] * qy[3] + qy[1] * qy[1];
    double invdet = 1 / det;
    qxz[0] = (qy[3] * q[0] - qy[1] * q[2]) * invdet;
    qxz[1] = 0;
    qxz[2] = (qy[1] * q[0] + qy[3] * q[2]) * invdet;
    //
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

double Get_Hip_Magnet(const double* values_pre, const double* values_cur)
{
    double x = values_cur[0] - values_pre[0];
    double y = values_cur[1] - values_pre[1];
    double z = values_cur[2] - values_pre[2];
    return sqrt(x * x + y * y + z * z);
}

void Gene_Hip_Match_Dataset(
    int nframe,
    std::vector<double>& hip_match_date,
    const std::vector<double>& timeseries_data,
    int channelSize)
{
    for (unsigned int i = 0; i < nframe; i++)
    {
        if (i == 0)
        {
            hip_match_date.push_back(0.0f);
        }
        else
        {
            double x = timeseries_data[i * channelSize + 0] - timeseries_data[(i - 1) * channelSize + 0];
            double y = timeseries_data[i * channelSize + 1] - timeseries_data[(i - 1) * channelSize + 1];
            double z = timeseries_data[i * channelSize + 2] - timeseries_data[(i - 1) * channelSize + 2];
            double speed_mag = sqrt(x * x + y * y + z * z);
            hip_match_date.push_back(speed_mag);
        }
    }
}


dfm2::CVec2d Catmul_spline(std::vector<dfm2::CVec2d>& p,float t)
{
    int p1, p2, p3, p0;
    p1 = (int)t + 1;
    p2 = p1 + 1;
    p3 = p2 + 1;
    p0 = p1 - 0;

    t = t - (int)t;

    float tt = t * t;
    float ttt = tt * t;

    float q1 = -ttt + 2.0f * tt - t;
    float q2 = 3.0f * ttt - 5.0f * tt + 2.0f;
    float q3 = -3.0f * ttt + 4.0f * tt + t;
    float q4 = ttt - tt;

    float tx = 0.5f * (p[p0].x * q1 + p[p1].x * q2 + p[p2].x * q3 + p[p3].x * q4);
    float ty = 0.5f * (p[p0].y * q1 + p[p1].y * q2 + p[p2].y * q3 + p[p3].y * q4);

    dfm2::CVec2d out{ tx,ty };
    return out;
}

void draw_spline(std::vector<dfm2::CVec2d>& p)
{
    ::glLineWidth(0.5);
    ::glColor3d(0.1, 0.1, 0.1);
    ::glBegin(GL_LINE_STRIP);
    for (auto k : p) 
    {
        ::glVertex3d(k.x, 0.1, k.y);
    }
    /*
    for (float t = 0.0; t < (float)p.size() - 3.0f; t += 0.01) {
        dfm2::CVec2d out = Catmul_spline(p, t);
        ::glVertex3d(out.x, 0.1, out.y);
    }
    */
    ::glEnd();
}



dfm2::CVec2d beizer_curve(dfm2::CVec2d& p1, dfm2::CVec2d& p2, dfm2::CVec2d& p3, double t)
{
    dfm2::CVec2d out{ (1 - t) * (1 - t) * p1.x + 2 * (1 - t) * t * p2.x + t * t * p3.x, (1 - t) * (1 - t) * p1.y + 2 * (1 - t) * t * p2.y + t * t * p3.y };
    return out;
}

void draw_curve(dfm2::CVec2d& p1, dfm2::CVec2d& p2, dfm2::CVec2d& p3)
{

        ::glLineWidth(1);
        ::glColor3d(1.0, 0, 0);
        ::glBegin(GL_LINE_STRIP);
        for (double j = 0; j < 1; j += 0.1) {
            dfm2::CVec2d out = beizer_curve(p1, p2, p3, j);
            ::glVertex3d(out.x, 0.1, out.y);
        }
        ::glEnd();

}

void draw_motion_line(std::vector<double>& motion_line)
{
    if (motion_line.size() < 18)
    {

    }
    else
    {
        for (int j = 0; j < 9; j += 3) {
            ::glLineWidth(1);
            ::glColor3d(1.0, 0, 0);
            ::glBegin(GL_LINE_STRIP);
            for (unsigned int i = j; i < motion_line.size(); i += 9) {
                ::glVertex3d(motion_line[i], motion_line[i + 1], motion_line[i + 2]);
            }
            ::glEnd();
        }
        
    }
}

void SetPose_BioVisionHierarchy_Rotate(
    std::vector<dfm2::CRigBone>& bones,
    const std::vector<dfm2::CChannel_BioVisionHierarchy>& channels,
    const double* values,
    dfm2::CVec2d& traj,
    dfm2::CVec2d& face_dir) {
    for (auto& bone : bones) {
        bone.quatRelativeRot[0] = 0.0;
        bone.quatRelativeRot[1] = 0.0;
        bone.quatRelativeRot[2] = 0.0;
        bone.quatRelativeRot[3] = 1.0;
    }
    const size_t nch = channels.size();
    int traj_index = 0;
    for (unsigned int ich = 0; ich < nch; ++ich) {
        const int ibone = channels[ich].ibone;
        const int iaxis = channels[ich].iaxis;
        const bool isrot = channels[ich].isrot;
        const double val = values[ich];
        assert(ibone < (int)bones.size());
        assert(iaxis >= 0 && iaxis < 3);
        if (!isrot) {
            if (traj_index == 0)
                bones[ibone].transRelative[iaxis] = traj.x;
            else if (traj_index == 2)
                bones[ibone].transRelative[iaxis] = traj.y;
            else
                bones[ibone].transRelative[iaxis] = val;
            traj_index += 1;
        }
        else if ((ich == 3) || (ich == 4) || (ich == 5))
        {
            if (ich == 3)
            {
                dfm2::CMat3d r = dfm2::CMat3d::Identity();
                for (unsigned int idim = 3; idim < 6; ++idim) {
                    const double ar = val * M_PI / 180.0;
                    double a[3] = { 0, 0, 0 };
                    a[channels[idim].iaxis] = ar;
                    dfm2::CMat3d dr;
                    dr.SetRotMatrix_Cartesian(a);
                    r = r * dr; // this order is important
                }
                dfm2::CQuatd q;
                //r.GetQuat_RotMatrix(q.p);
                q = r.GetQuaternion();
                dfm2::CQuatd qy, qxy;
                SeparateYRot(qy.p, qxy.p, q.p);
                dfm2::CQuatd qyinv = qy.conjugate();
                dfm2::CQuatd unroy_q;
                dfm2::QuatQuat(unroy_q.p, qyinv.p, q.p);
                double qtmp[4];
                dfm2::QuatQuat(qtmp,
                    bones[ibone].quatRelativeRot, unroy_q.p);


                dfm2::CQuatd q0y;
                const dfm2::CMat3d R0(
                    face_dir.y, 0, face_dir.x,
                    0, 1, 0,
                    -face_dir.x, 0, face_dir.y);
                q0y = R0.GetQuaternion();
                q0y.SetSmallerRotation();
                //dfm2::CQuatd q0yinv = q0y.conjugate();
                double qtmp2[4];
                dfm2::QuatQuat(qtmp2, q0y.p, qtmp);

                dfm2::Copy_Quat(bones[ibone].quatRelativeRot, qtmp2);
            }
        }
        else {
            const double ar = val * M_PI / 180.0;
            double v0[3] = { 0, 0, 0 };
            v0[iaxis] = 1.0;
            double dq[4] = { v0[0] * sin(ar * 0.5), v0[1] * sin(ar * 0.5), v0[2] * sin(ar * 0.5), cos(ar * 0.5) };
            double qtmp[4];
            dfm2::QuatQuat(qtmp,
                bones[ibone].quatRelativeRot, dq);
            /*
            if (ich == 5)
            {
                double final_in[4];
                double final_in_dir[4];
                dfm2::CQuatd qy, qxz;
                SeparateYRot(qy.p, qxz.p, qtmp);
                dfm2::CQuatd qyinv = qy.conjugate();
                dfm2::QuatQuat(final_in, qtmp, qyinv.p);
                dfm2::CQuatd q0y;
                const dfm2::CMat3d R0(
                    dirz.y, 0, dirz.x,
                    0, 1, 0,
                    -dirz.x, 0, dirz.y);
                q0y = R0.GetQuaternion();
                q0y.SetSmallerRotation();
                dfm2::CQuatd q0yinv = q0y.conjugate();
                dfm2::QuatQuat(final_in_dir, final_in, q0yinv.p);
                dfm2::Copy_Quat(bones[ibone].quatRelativeRot, final_in_dir);
            }
            */

            dfm2::Copy_Quat(bones[ibone].quatRelativeRot, qtmp);
          
        }
    }
    UpdateBoneRotTrans(bones);
}


int Best_match(std::vector<dfm2::CRigBone>& abone, 
    std::vector<dfm2::CChannel_BioVisionHierarchy> aChannel,
    std::vector<double> vec_bvh_time_series_data,
    std::vector<int> icandidate,
    dfm2::CVec2d root_pos)
{
    double cost = 600.0;
    int index_target = -1;
    int nch = aChannel.size();
    std::vector<double> curr_feature;
    for (unsigned int i = 0; i < abone.size(); i++)
    {
        auto bpos = abone[i].RootPosition();
        curr_feature.push_back(bpos[0]);
        curr_feature.push_back(bpos[1]);
        curr_feature.push_back(bpos[2]);
    }
    for (int ic = 0; ic < icandidate.size(); ic++)
    {
        double tempt_cost = 0.0f;
       // SetPose_BioVisionHierarchy_Rotate(abone, aChannel, vec_bvh_time_series_data.data() + icandidate[ic] * nch, root_pos);
        std::vector<double> pred_feature;
        for (unsigned int i = 0; i < abone.size(); i++)
        {
            auto bpos = abone[i].RootPosition();
            pred_feature.push_back(bpos[0]);
            pred_feature.push_back(bpos[1]);
            pred_feature.push_back(bpos[2]);
        }
        for (unsigned int k = 0; k < curr_feature.size(); k++)
        {
            tempt_cost += (pred_feature[k] - curr_feature[k]) * (pred_feature[k] - curr_feature[k]);
        }
        //std::cout << tempt_cost << std::endl;
        if (tempt_cost < cost)
        {
            index_target = ic;
            cost = tempt_cost;
        }
    }
    if (index_target == -1)
        return -1;
    else
        return icandidate[index_target];
}


void Gene_feature_vector(std::vector<dfm2::CRigBone>& pre_bone, std::vector<dfm2::CRigBone>& curr_bone,std::vector<double>& feature_vec)
{
    for (unsigned int i = 0; i < curr_bone.size(); i++)
    {
        auto pre_pos = pre_bone[i].RootPosition();
        auto pos = curr_bone[i].RootPosition();
        feature_vec[i * 3 + 0] = pos[0] - pre_pos[0];
        feature_vec[i * 3 + 1] = pos[1] - pre_pos[1];
        feature_vec[i * 3 + 2] = pos[2] - pre_pos[2];
    }
}


int Pose_Match_Best(std::vector<double>& pose_match_data, 
    std::vector<int>& candidate_index, 
    std::vector<double>& current_feature,
    int b_size)
{
    int best_index = 0;
    int f_size = b_size * 3;
    double cost = 9999999.0;
    for (unsigned int i = 0; i < candidate_index.size(); i++)
    {
        double tempt_cost = 0.0;
        int m_index = candidate_index[i] * f_size;
        for (unsigned int j = 0; j < f_size; j++)
        {
            tempt_cost += abs(current_feature[j] - pose_match_data[m_index + j]);
        }
        if (tempt_cost < cost)
        {
            best_index = i;
            cost = tempt_cost;
        }
    }
    return best_index;
}


void Gene_Pose_Match_Dataset(std::vector<dfm2::CRigBone>& abones, 
    std::vector<dfm2::CChannel_BioVisionHierarchy>& achannels, 
    int nframe,
    std::vector<double>& pose_match_date, 
    const std::vector<double>& timeseries_date, 
    int channelSize)
{
    int b_size = abones.size();
    std::vector<double> pre_abone(b_size * 3,0.0);
    for (unsigned int i = 0; i < nframe; i++)
    {
        SetPose_BioVisionHierarchy(abones, achannels, timeseries_date.data() + i * channelSize);
        for (unsigned int j = 0; j < b_size; j++)
        {
            auto pos = abones[j].RootPosition();
            if ((i == 0) && (j == 0))
            {
                pose_match_date.push_back(0.0f);
                pose_match_date.push_back(0.0f);
                pose_match_date.push_back(0.0f);
                pre_abone[j * 3 + 0] = pos[0];
                pre_abone[j * 3 + 1] = pos[1];
                pre_abone[j * 3 + 2] = pos[2];
            }
            else {
                double x_diff = pos[0] - pre_abone[j * 3 + 0];
                double y_diff = pos[1] - pre_abone[j * 3 + 1];
                double z_diff = pos[2] - pre_abone[j * 3 + 2];
                pose_match_date.push_back(x_diff);
                pose_match_date.push_back(y_diff);
                pose_match_date.push_back(z_diff);
                // std::cout << j << std::endl;
                 //std::cout << x_diff << "]" << y_diff << "]" << z_diff << std::endl;
                pre_abone[j * 3 + 0] = pos[0];
                pre_abone[j * 3 + 1] = pos[1];
                pre_abone[j * 3 + 2] = pos[2];
            }
        }
    }
}


std::vector<int> Traj_Match_Best_hip(std::vector<double>& traj_match_date,
    std::vector<double>& desired_future_traj,
    int database_size,
    int num_sample,
    std::vector<double>& hip_match_data,
    double target_speed_mag,
    int pre_frame)
{
    int low_bind = pre_frame - 120;
    float loss_fra = 20.f;
    int candid_size = 5;
    std::vector<double> cost(candid_size, 9999999.0);
    std::vector<int> target_frame(candid_size, 0);
    int total_frame = database_size / num_sample;
    //std::cout << "total" << total_frame << std::endl;
    for (unsigned int iframe = 0; iframe < total_frame; iframe++)
    {
        double hip_cost = 0.0f;
        double loss = abs(hip_match_data[iframe] - target_speed_mag);
        if ((iframe < pre_frame) && (iframe > low_bind))
        {
            hip_cost += 999999.0;
        }
        else
        {
            hip_cost += loss;
        }
        double tempt_cost = loss_fra * hip_cost;
        for (unsigned int isample = 0; isample < num_sample; isample += 2)
        {

            double x_diff = traj_match_date[iframe * num_sample + isample + 0] - desired_future_traj[isample + 0];
            double z_diff = traj_match_date[iframe * num_sample + isample + 1] - desired_future_traj[isample + 1];
            tempt_cost += sqrt(x_diff * x_diff + z_diff * z_diff);
        }
        for (unsigned int z = 0; z < candid_size; z++)
        {
            if (tempt_cost < cost[z])
            {
                target_frame[z] = iframe;
                cost[z] = tempt_cost;
                break;
            }
        }
    }
    return target_frame;
}


std::vector<int> Traj_Match_Best(std::vector<double>& traj_match_date,
    std::vector<double>& desired_future_traj,
    int database_size,
    int num_sample,
    int& c_f,
    double& ongoing_loss,
    int cur_frame)
{
    float threshold = 100;
    float loss_fra = 20.f;
    int candid_size = 1;
    std::vector<double> cost(candid_size,9999999.0);
    std::vector<int> target_frame(candid_size, 0);
    int total_frame = database_size / num_sample;
    //std::cout << "total" << total_frame << std::endl;
    for (unsigned int iframe = 0; iframe < total_frame; iframe++)
    {
        double tempt_cost = 0.0;
        for (unsigned int isample = 0; isample < num_sample; isample += 2)
        {

            double x_diff = traj_match_date[iframe * num_sample + isample + 0] - desired_future_traj[isample + 0];
            double z_diff = traj_match_date[iframe * num_sample + isample + 1] - desired_future_traj[isample + 1];
            tempt_cost += sqrt(x_diff * x_diff + z_diff * z_diff);
        }
        for (unsigned int z = 0; z < candid_size; z++)
        {
            if (tempt_cost < cost[z])
            {
                target_frame[z] = iframe;
                cost[z] = tempt_cost;
                break;
            }
        }
    }

    //return target_frame;
    
    if ((cost[0] + threshold) > ongoing_loss)
    {
        if (c_f == 100)
        {
            c_f = 0;
            ongoing_loss = cost[0];
            return target_frame;
        }
        else
        {
            c_f += 1;
            target_frame[0] = cur_frame + 1;
            return target_frame;
        }
    }
    else
    {
        c_f = 0;
        ongoing_loss = cost[0];
        return target_frame;
    }
    
}


void Gene_Traj_Match_Dataset(std::vector<double>& traj_match_data,const std::vector<double>& timeseries_date, std::vector<dfm2::CChannel_BioVisionHierarchy>& aChannelRotTransBone, int sample_frame)
{
    const unsigned int channelSize = aChannelRotTransBone.size();
    int size_date = timeseries_date.size();
    int size_timeseries = size_date / channelSize;
    for (unsigned int ich = 0; ich < size_timeseries; ich += 1)
    {
        for (unsigned int j = 1; j < 6; j += 1)
        {
            //const double val = (*timeseries_date)[channelSize * ich + j * sample_frame + 0];
            double x_root = timeseries_date[channelSize * ich + 0];
            double z_root = timeseries_date[channelSize * ich + 2];
            int ix = channelSize * ich + j * sample_frame * channelSize + 0;
            int iz = channelSize * ich + j * sample_frame * channelSize + 2;

            dfm2::CMat3d r = dfm2::CMat3d::Identity();
            for (unsigned int idim = 0; idim < 3; ++idim) {
                double a[3] = { 0, 0, 0 };
                const unsigned int iaxis = aChannelRotTransBone[3 + idim].iaxis;
                a[iaxis] = timeseries_date[ich * channelSize + 3 + idim] * M_PI / 180.0;
                dfm2::CMat3d dr;
                dr.SetRotMatrix_Cartesian(a);
                r = r * dr; // this order is important
            }
            dfm2::CQuatd q;
            q = r.GetQuaternion();
            dfm2::CQuatd qy, qxz;
            SeparateYRot(qy.p, qxz.p, q.p);
            dfm2::CQuatd qyinv = qy.conjugate();

            if (ix >= size_date)
            {
                traj_match_data.push_back(0.0f);
                traj_match_data.push_back(0.0f);
            }
            else 
            {
                double x_pos = timeseries_date[ix] - x_root;
                double y_pos = timeseries_date[iz] - z_root;
                double traj_local[3] = { x_pos , 0.0, y_pos };
                double un_Y_traj[3];
                dfm2::QuatVec(un_Y_traj, qyinv.p, traj_local);
                traj_match_data.push_back(un_Y_traj[0]);
                traj_match_data.push_back(un_Y_traj[2]);
            }
        }
    }
}


void Read_BioVisionHierarchy_skip(
    std::vector<dfm2::CRigBone>& bones,
    std::vector<dfm2::CChannel_BioVisionHierarchy>& channels,
    size_t& nframe,
    double& frame_time,
    std::vector<double>& frame_channel,
    std::string& bvh_header,
    const std::string& path_bvh,
    int skip) {
    std::ifstream fin;
    fin.open(path_bvh.c_str());
    if (!fin.is_open()) {
        std::cout << "cannot open file" << std::endl;
        return;
    }
    bvh_header.clear();
    bones.clear();
    channels.clear();
    //
    std::string line;
    std::vector<int> stackIndBone;
    while (std::getline(fin, line)) {
        if (line[line.size() - 1] == '\n') line.erase(line.size() - 1); // remove the newline code
        if (line[line.size() - 1] == '\r') line.erase(line.size() - 1); // remove the newline code
        bvh_header += (line + '\n');
        line = dfm2::rig_bvh::MyReplace(line, '\t', ' ');
        std::vector<std::string> aToken = dfm2::rig_bvh::MySplit(line, ' ');
        //    std::cout << aToken[0] << std::endl;
        if (aToken[0] == "HIERARCHY") {
            assert(bones.empty());
        }
        else if (aToken[0] == "ROOT") {
            assert(bones.empty());
            dfm2::CRigBone br;
            assert(aToken.size() == 2);
            br.name = aToken[1];
            bones.push_back(br);
        }
        else if (aToken[0] == "{") {
            stackIndBone.push_back(static_cast<int>(bones.size() - 1));
            if (stackIndBone.size() > 1) {
                int ibp = stackIndBone[stackIndBone.size() - 2];
                auto ib = static_cast<unsigned int>(bones.size() - 1);
                bones[ib].ibone_parent = ibp;
            }
        }
        else if (aToken[0] == "}") {
            stackIndBone.resize(stackIndBone.size() - 1);
        }
        else if (aToken[0] == "OFFSET") {
            assert(aToken.size() == 4);
            size_t ib = bones.size() - 1;
            double org_x = dfm2::rig_bvh::myStod(aToken[1]);
            double org_y = dfm2::rig_bvh::myStod(aToken[2]);
            double org_z = dfm2::rig_bvh::myStod(aToken[3]);
            bones[ib].invBindMat[3] = -org_x;
            bones[ib].invBindMat[7] = -org_y;
            bones[ib].invBindMat[11] = -org_z;
            if (stackIndBone.size() > 1) {
                const int ibp = stackIndBone[stackIndBone.size() - 2];
                assert(ibp < (int)bones.size());
                bones[ib].invBindMat[3] += bones[ibp].invBindMat[3];
                bones[ib].invBindMat[7] += bones[ibp].invBindMat[7];
                bones[ib].invBindMat[11] += bones[ibp].invBindMat[11];
            }
        }
        else if (aToken[0] == "CHANNELS") {
            assert(aToken.size() >= 2);
            int nch = dfm2::rig_bvh::myStoi(aToken[1]);
            assert((int)aToken.size() == nch + 2);
            assert(!bones.empty());
            const auto ib = static_cast<unsigned int>(bones.size() - 1);
            for (int ich = 0; ich < nch; ++ich) {
                const std::string& type_ch = aToken[ich + 2];
                if (type_ch == "Xposition") { channels.emplace_back(ib, 0, false); }
                else if (type_ch == "Yposition") { channels.emplace_back(ib, 1, false); }
                else if (type_ch == "Zposition") { channels.emplace_back(ib, 2, false); }
                else if (type_ch == "Xrotation") { channels.emplace_back(ib, 0, true); }
                else if (type_ch == "Yrotation") { channels.emplace_back(ib, 1, true); }
                else if (type_ch == "Zrotation") { channels.emplace_back(ib, 2, true); }
                else {
                    std::cout << "ERROR-->undefiend type" << std::endl;
                }
            }
        }
        else if (aToken[0] == "JOINT") {
            dfm2::CRigBone br;
            assert(aToken.size() == 2);
            br.name = aToken[1];
            bones.push_back(br);
        }
        else if (aToken[0] == "End") {
            assert(aToken[1] == "Site");
            dfm2::CRigBone br;
            assert(aToken.size() == 2);
            br.name = aToken[1];
            bones.push_back(br);
        }
        else if (aToken[0] == "MOTION") {
            break;
        }
    }
    nframe = 0;
    {
        std::string stmp0;
        std::getline(fin, line);  // Frames: ***
        std::stringstream ss(line);
        ss >> stmp0 >> nframe;
        // std::cout << "frame: " << nframe << std::endl;
    }
    {
        std::string stmp0, stmp1;
        std::getline(fin, line);  // Frame Time: ***
        std::stringstream ss(line);
        ss >> stmp0 >> stmp1 >> frame_time;
        // std::cout << "frametime: " << frame_time << std::endl;
    }
    nframe = nframe - skip;
    const size_t nchannel = channels.size();
    frame_channel.resize(nframe * nchannel);
    for (unsigned int i = 0; i < skip; i++)
    {
        std::getline(fin, line);
    }
    for (unsigned int iframe = 0; iframe < nframe; ++iframe) {
        std::getline(fin, line);
        line = dfm2::rig_bvh::MyReplace(line, '\t', ' ');
        if (line[line.size() - 1] == '\n') line.erase(line.size() - 1); // remove the newline code
        if (line[line.size() - 1] == '\r') line.erase(line.size() - 1); // remove the newline code
        std::vector<std::string> aToken = dfm2::rig_bvh::MySplit(line, ' ');
        //    std::cout << aToken.size() << " " << aChannelRotTransBone.size() << std::endl;
        assert(aToken.size() == channels.size());
        for (unsigned int ich = 0; ich < nchannel; ++ich) {
            frame_channel[iframe * nchannel + ich] = dfm2::rig_bvh::myStod(aToken[ich]);
        }
    }
    // ---------------
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


float x_target_v = 0.f;
float y_target_v = 0.f;
float keyboard_speed_mag = 25.0f;

int normal_frame = 0;
bool set_y_0 = true;
int main(int argc, char* argv[]) {
    const int match_gap_frame = 2;
    const int sample_frame = 20;

    //const int boneSize= 38;
    /*
    std::vector<dfm2::CRigBone> readBone;
    readBone.resize(38);
    std::ifstream inFile;
    inFile.open(std::string(PATH_SOURCE_DIR) + "/../data/motion_match_data_1.bin", std::ios::in | std::ios::binary);

    inFile.read((char*)&readBone, readBone.size());
    if (inFile.fail()) { std::cout << "Load Binary Error" << std::endl; return false; }
    std::cout << "size of current bone" << std::endl;
    std::cout << readBone.size() << std::endl;
    */



    std::vector<std::string> vec_path_bvh = { 
        "/../data/LocomotionFlat01_000.bvh",
        "/../data/LocomotionFlat01_000_mirror.bvh",
        "/../data/LocomotionFlat02_000.bvh",
         "/../data/LocomotionFlat02_000_mirror.bvh",
        "/../data/LocomotionFlat03_000.bvh",
         "/../data/LocomotionFlat03_000_mirror.bvh",
        "/../data/LocomotionFlat04_000.bvh",
         "/../data/LocomotionFlat04_000_mirror.bvh",
        "/../data/LocomotionFlat05_000.bvh",
         "/../data/LocomotionFlat05_000_mirror.bvh",
        "/../data/LocomotionFlat06_000.bvh",
         "/../data/LocomotionFlat06_000_mirror.bvh",
        "/../data/LocomotionFlat07_000.bvh"
         "/../data/LocomotionFlat07_000_mirror.bvh"
    };

    std::vector<dfm2::CRigBone> aBone;
    std::vector<dfm2::CChannel_BioVisionHierarchy> aChannelRotTransBone;
    double frame_time;
    std::string header_bvh;
    std::vector<double> traj_match_data;
    std::vector<double> pose_match_data;
    std::vector<double> vec_bvh_time_series_data;
    std::vector<double> hip_match_date;
    
    for (auto& ipath : vec_path_bvh) {
        std::vector<double> series_data_tempt;
        size_t nframe = 0;
        std::string path_bvh = std::string(PATH_SOURCE_DIR) + ipath;
        Read_BioVisionHierarchy_skip(
            aBone,
            aChannelRotTransBone,
            nframe,
            frame_time,
            series_data_tempt,
            header_bvh,
            path_bvh, 50);

        std::vector<double> traj_match_data_tempt;
        std::vector<double> pose_match_data_tempt;
        std::vector<double> hip_match_date_tempt;
        
        Gene_Traj_Match_Dataset(traj_match_data_tempt, series_data_tempt, aChannelRotTransBone, sample_frame);
        //Gene_Pose_Match_Dataset(aBone, aChannelRotTransBone, nframe, pose_match_data_tempt, series_data_tempt, 96);

        vec_bvh_time_series_data.insert(vec_bvh_time_series_data.end(), series_data_tempt.begin(), series_data_tempt.end());
        Gene_Hip_Match_Dataset(nframe, hip_match_date_tempt, vec_bvh_time_series_data, 96);
        traj_match_data.insert(traj_match_data.end(), traj_match_data_tempt.begin(), traj_match_data_tempt.end());
        //pose_match_data.insert(pose_match_data.end(), pose_match_data_tempt.begin(), pose_match_data_tempt.end());

        hip_match_date.insert(hip_match_date.end(), hip_match_date_tempt.begin(), hip_match_date_tempt.end());
    }
        /*
    if (set_y_0 == true)
    {
        for (unsigned int i = 0; i < vec_bvh_time_series_data.size() / 96; i++)
        {
            double thex = abs(vec_bvh_time_series_data[i * 96 + 3]);
            double they = vec_bvh_time_series_data[i * 96 + 4];
            if (thex > 160)
            {
                vec_bvh_time_series_data[i * 96 + 4] = 180.0;
            }
            else
            {
                vec_bvh_time_series_data[i * 96 + 4] = 0.0;
            }
            
            vec_bvh_time_series_data[i * 96 + 0] = 0.0;
            vec_bvh_time_series_data[i * 96 + 2] = 0.0;
        }
    }
    */
        dfm2::UpdateBoneRotTrans(aBone);
        std::cout << "aBone: " << aBone.size() << std::endl;
        std::cout << "aChannels: " << aBone.size() << std::endl;
        assert(aBone.size() == 38 && aChannelRotTransBone.size() == 96);
        std::cout << "traj_match_data: " << traj_match_data.size() / 10 << std::endl;
        std::cout << "vec_bvh_time_series_data: " << vec_bvh_time_series_data.size()/96 << std::endl;
        std::cout << "pose_match_data: " << pose_match_data.size() / (3 * 38) << std::endl;
        const int total_f = vec_bvh_time_series_data.size() / 96;

        // 
        // 
        //while (!fin.eof()) {
        //    double a = -1;
        //    fin >> a;
        //    if (!fin) { break; }  // the line break at the end of the file
        //    vec_phase.push_back(a);
        //}
        /*
        {
            std::ofstream outFile;
            std::string name = "motion_match_data_1";
            outFile.open(
                std::string(PATH_SOURCE_DIR) +
                "/../data/" + name +
                ".bin", std::ios::binary);
            outFile.write(
                reinterpret_cast<const char*>(aBone.data()),
                sizeof(double) * aBone.size());
            outFile.close();
        }
        */
        //exit(3);
    class MytempViewer : public delfem2::glfw::CViewer3 {
    public:
        MytempViewer() : CViewer3(45) {
        }

        void key_press(int key, int mods) override {
            delfem2::glfw::CViewer3::key_press(key, mods);
            if (key == GLFW_KEY_F) {
                if (keyboard_speed_mag == 25.0f)
                    keyboard_speed_mag = 8.0f;
                else
                    keyboard_speed_mag = 25.0f;
            }
            if (key == GLFW_KEY_X) {
                x_target_v = 0.0f;
                y_target_v = 0.f;
            }
        }

        void key_repeat(int key, int mods) override {
            delfem2::glfw::CViewer3::key_press(key, mods);
            if (key == GLFW_KEY_W) {
                x_target_v = 0.f;
                y_target_v = -keyboard_speed_mag;
            }
            if (key == GLFW_KEY_S) {
                x_target_v = 0.f;
                y_target_v = keyboard_speed_mag;
            }
            if (key == GLFW_KEY_A) {
                x_target_v = -keyboard_speed_mag;
                y_target_v = 0.f;
            }
            if (key == GLFW_KEY_D) {
                x_target_v = keyboard_speed_mag;
                y_target_v = 0.f;
            }
        }
    };

    MytempViewer viewer_source;

    viewer_source.width = 1920;
    viewer_source.height = 1080;
    viewer_source.window_title = "Motion Matching";
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

    int path_count = -1;

    static float f = 0.01f;
    double speed_x = 0.0f;
    double speed_y = 0.0f;
    //damper_p the_p = damper_p();

    enum
    {
        TRAJ_MAX = 120,
        TRAJ_SUB = 20,
        PRED_MAX = 6,
        PRED_SUB = 20,
    };

    float trajx_prev[TRAJ_MAX];
    float trajy_prev[TRAJ_MAX];

    float predx[PRED_MAX], predy[PRED_MAX];
    float predxv[PRED_MAX], predyv[PRED_MAX];
    float predxa[PRED_MAX], predya[PRED_MAX];

    float halflife = 0.70f;
    float dt = 1.0 / 60.0f;
    //float timescale = 240.0f;
    float trajx = 0.0f;
    float trajy = 0.0f;
    float trajxv = 0.0, trajyv = 0.0;
    float trajxa = 0.0, trajya = 0.0;
    float traj_xv_goal = 0.0;
    float traj_yv_goal = 0.0;


    for (int i = 0; i < TRAJ_MAX; i++)
    {
        trajx_prev[i] = 0.0f;
        trajy_prev[i] = 0.0f;
    }

    float velocity_mag = 10.0f;
    dfm2::CVec2d face_dirZ(1.0f, 0.f);
    static int iframe = 0;
    const int nch = aChannelRotTransBone.size();
    std::vector<double> current_feature(38 * 3, 0.0);
    std::vector<dfm2::CRigBone> pre_bone;
    pre_bone = aBone;

    int frame_check = 0;
    int current_ani_index = 0;

    std::vector<double> motion_line;

    int pre_f = 0;
    int curr_f = 0;

    int c_f = 0;
    double ongoing_loss = 9999999.0;

    while (!glfwWindowShouldClose(viewer_source.window))
    {
        glfwMakeContextCurrent(viewer_source.window);

        if (frame_check == match_gap_frame)
        {
            frame_check = 0;
        }
        bool facing_control_mode = false;


        Floor floor{ 100, +0.1 };


        viewer_source.DrawBegin_oldGL();
        ::glEnable(GL_LINE_SMOOTH);
        ::glEnable(GL_BLEND);
        ::glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);



        for (int i = TRAJ_MAX - 1; i > 0; i--)
        {
            trajx_prev[i] = trajx_prev[i - 1];
            trajy_prev[i] = trajy_prev[i - 1];
        }



        float gamepad_x = 0.0f;
        float gamepad_y = 0.0f;
        float gamepad2_x = 0.0f;
        float gamepad2_y = 0.0f;

        bool buttonX = false;
        bool buttonY = false;
        bool buttonA = false;
        bool buttonB = false;


        // input bind with controller ++++++++++++
        bool check_controller = false;
        int present = 0;
        
        double x_mag = 0.0f;
        double y_mag = 0.0f;
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
                //std::cout << "Right X : " << axes[2] << std::endl;
                //std::cout << "Right Y : " << axes[3] << std::endl;
                //std::cout << "trigger LT : " << axes[4] << std::endl;
                //std::cout << "trigger RT : " << axes[5] << std::endl;
                //std::cout << "trigger R3 : " << axes[6] << std::endl;
                //std::cout << "trigger R4 : " << axes[7] << std::endl;

                if (gamepad_mag > 0.2f)
                {
                    float gamepaddirx = gamepad_x / gamepad_mag;
                    float gamepaddiry = gamepad_y / gamepad_mag;
                    float gamepadclippedmag = gamepad_mag > 1.0f ? 1.0f : gamepad_mag * gamepad_mag;
                    gamepad_x = gamepaddirx * gamepadclippedmag;
                    gamepad_y = gamepaddiry * gamepadclippedmag;

                }
                else
                {
                    gamepad_x = 0.0f;
                    gamepad_y = 0.0f;
                }
                //22.0f close to run speed 
                //
                if (triggerLT == 1) {
                    facing_control_mode = true;
                }
                if (triggerRT == 1) {
                    // run

                    f = damper(f, 0.01f, 0.1f);
                    //f = 0.01f;
                    halflife = damper(halflife, 1.0f, 0.1f);
                    //halflife = 0.7f;
                    velocity_mag = damper(velocity_mag, 25.0f, 0.1f);
                    traj_xv_goal = gamepad_x * velocity_mag;
                    traj_yv_goal = gamepad_y * velocity_mag;
                }
                else
                {
                    // walk 

                    f = damper(f, 0.006f, 0.1f);
                    //f = 0.006f;
                    halflife = damper(halflife, 0.9f, 0.1f);
                    //halflife = 0.4f;
                    velocity_mag = damper(velocity_mag, 10.0f, 0.1f);
                    traj_xv_goal = gamepad_x * velocity_mag;
                    traj_yv_goal = gamepad_y * velocity_mag;
                }
            }

            int buttonCont;
            const unsigned char* buttons = glfwGetJoystickButtons(GLFW_JOYSTICK_1, &buttonCont);
            //0 = A, 1 = B, 2 = X, 3 = Y
            if (GLFW_PRESS == buttons[0])
            {
                buttonA = true;
            }
            else if (GLFW_RELEASE == buttons[0])
            {
                buttonA = false;
            }
            if (GLFW_PRESS == buttons[1])
            {
                buttonB = true;
            }
            else if (GLFW_RELEASE == buttons[1])
            {
                buttonB = false;
            }
            if (GLFW_PRESS == buttons[2])
            {
                buttonX = true;
            }
            else if (GLFW_RELEASE == buttons[2])
            {
                buttonX = false;
            }
            if (GLFW_PRESS == buttons[3])
            {
                buttonY = true;
            }
            else if (GLFW_RELEASE == buttons[3])
            {
                buttonY = false;
            }
        }
        else
        {
            traj_xv_goal = damper(traj_xv_goal, x_target_v, 0.1f);
            traj_yv_goal = damper(traj_yv_goal, y_target_v, 0.1f);
        }


        spring_character_update(trajx, trajxv, trajxa, traj_xv_goal, halflife, dt);
        spring_character_update(trajy, trajyv, trajya, traj_yv_goal, halflife, dt);

        spring_character_predict(predx, predxv, predxa, PRED_MAX, trajx, trajxv, trajxa, traj_xv_goal, halflife, dt * PRED_SUB);
        spring_character_predict(predy, predyv, predya, PRED_MAX, trajy, trajyv, trajya, traj_yv_goal, halflife, dt * PRED_SUB);

        trajx_prev[0] = trajx;
        trajy_prev[0] = trajy;

        std::vector<dfm2::CVec2d> vec_pos2;
        vec_pos2.clear();
        std::vector<double> future_traj;
        future_traj.clear();

        

        for (int i = 0; i < TRAJ_MAX - TRAJ_SUB; i += TRAJ_SUB)
        {
            std::vector start = { trajx_prev[i + 0], trajy_prev[i + 0] };
            std::vector stop = { trajx_prev[i + TRAJ_SUB], trajy_prev[i + TRAJ_SUB] };

            if (i == 0) {
                draw_black_sphere(start);
                dfm2::CVec2d tempt(start[0], start[1]);
                vec_pos2.push_back(tempt);
                //std::cout << "root pos : " << start[0] << " - " << start[1] << std::endl;
            }
            else
            {
                dfm2::CVec2d tempt(start[0], start[1]);
                vec_pos2.push_back(tempt);
                draw_blue_sphere(start);
            }
            //DrawLineV(start, stop, BLUE);

            if (i + TRAJ_SUB == TRAJ_MAX - TRAJ_SUB)
            {
                dfm2::CVec2d tempt(stop[0], stop[1]);
                vec_pos2.push_back(tempt);
                draw_blue_sphere(stop);
            }
        }
        std::reverse(vec_pos2.begin(), vec_pos2.end());
    
        for (int i = 1; i < PRED_MAX; i += 1)
        {
            std::vector start = { predx[i + 0], predy[i + 0] };
            std::vector stop = { predx[i - 1], predy[i - 1] };

            //DrawLineV(start, stop, MAROON);
            dfm2::CVec2d tempt(start[0], start[1]);
            vec_pos2.push_back(tempt);
            //relative future traj
            future_traj.push_back(start[0] - trajx_prev[0]);
            future_traj.push_back(start[1] - trajy_prev[0]);
            draw_red_sphere(start);
        }


        speed_x = (predx[1] - trajx_prev[0]) / (dt * 20);
        speed_y = (predy[1] - trajy_prev[0]) / (dt * 20);
        if (facing_control_mode == false)
            viewer_source.view_rotation->Rot_Camera(-static_cast<float>(gamepad2_x), -static_cast<float>(gamepad2_y));


        //gamepad input are here  *** (ßÍß) ***
        float goal_x = 0.0f;
        float goal_y = 0.0f;
        dfm2::CVec2d goal_dirZ;
        if (facing_control_mode == false) {
            goal_dirZ[0] = predx[1] - trajx_prev[0];
            goal_dirZ[1] = predy[1] - trajy_prev[0];
            goal_dirZ.normalize();
            face_dirZ.normalize();
            goal_x = damper(face_dirZ.x, goal_dirZ.x, 0.1);
            goal_y = damper(face_dirZ.y, goal_dirZ.y, 0.1);
            if ((goal_x >= -1) && (goal_x <= 1))
            {
                face_dirZ.x = goal_x; // get rid of initial broken data, fix it later
            }
            if ((goal_y >= -1) && (goal_y <= 1))
            {
                face_dirZ.y = goal_y;
            }
        }
        else {
            goal_dirZ[0] = gamepad2_x;
            goal_dirZ[1] = gamepad2_y;
            goal_dirZ.normalize();
            face_dirZ.normalize();
            goal_x = damper(face_dirZ.x, goal_dirZ.x, 0.05);
            goal_y = damper(face_dirZ.y, goal_dirZ.y, 0.05);
            if ((goal_x >= -1) && (goal_x <= 1))
            {
                face_dirZ.x = goal_x; // get rid of initial broken data, fix it later
            }
            if ((goal_y >= -1) && (goal_y <= 1))
            {
                face_dirZ.y = goal_y;
            }
        }


        //face_dirZ.x = goal_x;
        //face_dirZ.y = goal_y;
        face_dirZ.normalize();
        //std::cout << "face_dir" << face_dirZ[0] << "-" << face_dirZ[1] << std::endl;
        //

        ::glLineWidth(1);
        ::glColor3d(0, 0, 0);
        ::glBegin(GL_LINES);
        ::glVertex3d(trajx_prev[0], 0.1f, trajy_prev[1]);
        ::glVertex3d(trajx_prev[0] + face_dirZ[0] * 5, 0.1f, trajy_prev[1] + face_dirZ[1] * 5);
        ::glEnd();

        dfm2::CVec2d root_pos2(trajx_prev[0], trajy_prev[0]);
        // std::cout <<   "traj size : " << vec_pos2.size() << std::endl;
                //  BVH test

        
        //Gene_feature_vector(pre_bone, aBone, current_feature);
        //int jump_frame = Pose_Match_Best(pose_match_data, jump_frame_candidate, current_feature, 38);

        //double hip_mag = Get_Hip_Magnet(vec_bvh_time_series_data.data() + pre_f * nch, vec_bvh_time_series_data.data() + curr_f * nch);

        //std::vector<int> jump_frame_candidate = Traj_Match_Best_hip(traj_match_data, future_traj, traj_match_data.size(), 10, hip_match_date, hip_mag, curr_f);
        std::vector<double> unY_traj;
        unrotate_Y_future_traj(future_traj, goal_dirZ, unY_traj);

        std::vector<int> jump_frame_candidate = Traj_Match_Best(traj_match_data, unY_traj, traj_match_data.size(), 10, c_f, ongoing_loss, curr_f);

        int jump_frame = jump_frame_candidate[0];
        curr_f = jump_frame;
        
        //if (jump_frame == curr_f)
        //    jump_frame += 1;
        //pre_f = curr_f;
        //curr_f = jump_frame;

        /*
        if (frame_check == 0)
        {
            std::vector<int> jump_frame_candidate = Traj_Match_Best(traj_match_data, future_traj, traj_match_data.size(), 10);
            int jump_frame = Best_match(aBone, aChannelRotTransBone, vec_bvh_time_series_data, jump_frame_candidate, root_pos2);
            if (jump_frame == -1)
            { 
                current_ani_index += 1;
                SetPose_BioVisionHierarchy_Rotate(aBone, aChannelRotTransBone, vec_bvh_time_series_data.data() + current_ani_index * nch, root_pos2);
            }
            else
            {
                current_ani_index = jump_frame;
                SetPose_BioVisionHierarchy_Rotate(aBone, aChannelRotTransBone, vec_bvh_time_series_data.data() + jump_frame * nch, root_pos2);
            }
        }
        else
        {
            if ((current_ani_index + 1) >= total_f)
            {
                SetPose_BioVisionHierarchy_Rotate(aBone, aChannelRotTransBone, vec_bvh_time_series_data.data() + current_ani_index * nch, root_pos2);
            }
            else 
            {
                current_ani_index += 1;
                SetPose_BioVisionHierarchy_Rotate(aBone, aChannelRotTransBone, vec_bvh_time_series_data.data() + current_ani_index * nch, root_pos2);
            }
        }
        */
        //trajx_prev[0] = vec_bvh_time_series_data[jump_frame * nch + 0];
        //trajy_prev[0] = vec_bvh_time_series_data[jump_frame * nch + 2];
        // 
        // 
        // ==========================DEBUG===========================
        //SetPose_BioVisionHierarchy(aBone, aChannelRotTransBone, vec_bvh_time_series_data.data() + normal_frame * nch);
        draw_coordinate();
        //debug_future_traj_correct(unY_traj,goal_dirZ);
        //debug_traj_correct(traj_match_data, normal_frame);
        //normal_frame += 1;
        // =========================DEBUG=========================
        // 
        // 
        SetPose_BioVisionHierarchy_Rotate(aBone, aChannelRotTransBone, vec_bvh_time_series_data.data() + jump_frame * nch, root_pos2,goal_dirZ);
        draw_spline(vec_pos2);
        //unY_traj.clear();
        /*
        {
            motion_line.push_back(aBone[0].RootPosition()[0]);
            motion_line.push_back(aBone[0].RootPosition()[1]);
            motion_line.push_back(aBone[0].RootPosition()[2]);
            motion_line.push_back(aBone[5].RootPosition()[0]);
            motion_line.push_back(aBone[5].RootPosition()[1]);
            motion_line.push_back(aBone[5].RootPosition()[2]);
            motion_line.push_back(aBone[21].RootPosition()[0]);
            motion_line.push_back(aBone[21].RootPosition()[1]);
            motion_line.push_back(aBone[21].RootPosition()[2]);
        }
        draw_motion_line(motion_line);
        */
        dfm2::opengl::DrawBone_Octahedron(
            aBone,
            5, -1,
            0.1, 1.0);
        
        // 
        DrawFloorShadow(aBone, floor, +0.1);

        floor.draw_checkerboard();
        frame_check++;
        // GUI *******
        ImGui_ImplOpenGL2_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        // 1. Show the big demo window (Most of the sample code is in ImGui::ShowDemoWindow()! You can browse its code to learn more about Dear ImGui!).
        {
            //static float f = 0.0f;
            ImGui::Begin("Guide");
            ImGui::Text("Controller User");
            ImGui::BulletText("Leftstick       - Control Moving Direction");
            ImGui::BulletText("Rightstick      - Control Camera Direction");
            ImGui::BulletText("LT + Rightstick - Control Facing Direction");
            ImGui::BulletText("RT (HOLD)       - Switch to Run");

            ImGui::Separator();
            ImGui::Text("Keyboard User");
            ImGui::BulletText(" W S A D (HOLD)                 - Control Moving Direction");
            ImGui::BulletText(" ALT (HOLD) + Left Click (HOLD) - Control Camera Direction");
            ImGui::BulletText(" F                              - Switch to Run/Walk");
            ImGui::BulletText(" X                              - Stop Character");
            ImGui::End();

            ImGui::Begin("Control Panel");                       
            ImGui::Checkbox("XBOX Gamepad", &check_controller);
            ImGui::Text("Buttons"); ImGui::SameLine();
            ImGui::Checkbox("X", &buttonX); ImGui::SameLine();
            ImGui::Checkbox("Y", &buttonY); ImGui::SameLine();
            ImGui::Checkbox("A", &buttonA); ImGui::SameLine();
            ImGui::Checkbox("B", &buttonB); ImGui::SameLine();
            ImGui::Separator();
            ImGui::Text("Trajectory smoothness");
            ImGui::SliderFloat("Traj Half-life", &halflife, 0.0f, 0.9f);
            ImGui::Separator();
            ImGui::Text("Frame_Life : %i", c_f);
            ImGui::Separator();
            ImGui::Text("Current Direction");
            ImGui::Text("dir x: %f ; dir y: %f", goal_x, goal_y);
            ImGui::Text("Current Speed");               // Display some text (you can use a format strings too)
            ImGui::Text("dir x: %f ; dir y: %f", speed_x, speed_y);
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