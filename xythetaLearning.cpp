#include <iostream>
#include <vector>
#include <cmath>
#include <math.h>
#include <cstdlib>
#include <algorithm>
#include <random>
#include <ctime>
#include <ode/ode.h>
#include <drawstuff/drawstuff.h>
#ifdef dDOUBLE                      
#define dsDrawSphere dsDrawSphereD 
#endif
#include "texturepath.h"
#define dsDrawBox dsDrawBoxD

#ifdef _MSC_VER
#pragma warning(disable:4244 4305)  // for VC++, no precision loss complaints
#endif

dWorldID world;  // 動力学計算用ワールド
dSpaceID space;  // 衝突検出用スペース
dGeomID  ground; // 地面
dGeomID slope;   // 斜面
dJointGroupID contactgroup; // コンタクトグループ
dReal stepsize = 0.01;  //ステップサイズ
dReal lx = 0.5,ly = 0.5, lz = 0.5, m  = 1.0;
dsFunctions fn;
static dReal side[3] = {lx, ly, lz};
FILE *fp, *fp3;
FILE *fp2;

typedef struct {       // MyObject構造体
  dBodyID body;        // ボディ(剛体)のID番号（動力学計算用）
  dGeomID geom;        // ジオメトリのID番号(衝突検出計算用）
  double lx, ly, lz, m;       // 縦，横，高さ，質量
} MyObject;
MyObject box[3];


// コールバック関数
static void nearCallback(void *data, dGeomID o1, dGeomID o2)
{
  static const int N = 10; // 接触点数の最大値
  dContact contact[N];     // 接触点

  int isGround = ((ground == o1) || (slope == o1))
               ||((ground == o2) || (slope == o2))
               ||((box[1].geom == o1) || (box[1].geom == o2))
               ||((box[0].geom == o1) || (box[0].geom == o2));

  // 衝突情報の生成 nは衝突点数
  int n =  dCollide(o1,o2,N,&contact[0].geom,sizeof(dContact));
  
     for (int i = 0; i < n; i++) {
       contact[i].surface.mode = dContactBounce; // 接触面の反発性を設定
       if (isGround) contact[i].surface.bounce = 0.02;          // 反発係数(0.0から1.0)
       else contact[i].surface.bounce = 0.0;
       contact[i].surface.bounce_vel = 0.0;      // 反発に必要な最低速度
       contact[i].surface.mu = 1.0;

       // 接触ジョイントの生成
       dJointID c = dJointCreateContact(world,contactgroup,
                                                &contact[i]);
       // 接触している２つの剛体を接触ジョイントにより拘束
       dJointAttach(c,dGeomGetBody(contact[i].geom.g1),
                      dGeomGetBody(contact[i].geom.g2));
    }
  
};

class Stickness{
private:
    double gravity;
    double massbox;
    double force_mag;
    double tau;
    double mu_static;
    double mu_dynamic;

    double x_threshold;
    double y_threshold;
    

    std::vector<double> state;
    const dReal *forceR;
    const dReal *forceL;
    int steps_beyond_terminated;
    
    const dReal *pos,*R,*R1,*R2;
    const dReal *force1,*force2,*force3,*T1,*T2,*T3,*Totalforce1,*Totalforce2;
    const dReal *pos1, *pos2,*pos3; 
    dReal relPos[3],relPos2[3];
    dReal fx_lim,fy_lim;

public:
    Stickness(){
        gravity = 9.8;
        massbox = 1.0;
        force_mag = 10.0;
        tau = 0.02;
        double mu_static = 0.8;
        double mu_dynamic = 0.5;

        x_threshold = 2.4;
        y_threshold = 2.4;

        state.resize(9);
        
        steps_beyond_terminated = 0;
    }

    bool terminate() {
        double x1 = state[0];
        double x2 = state[4];
        double y1 = state[2];
        double y2 = state[6];
        double xdistance = x2 -x1;
        double ydistance = y2 -y1;
        if (abs(xdistance) >2.0 || abs(ydistance) > 2.0){
            return true;
        }
        else{
            return false;
        }
    }
    void calculatestickness(dBodyID body1, dBodyID body2){
        fx_lim = 20;
        fy_lim = 20;
        
        pos1 = dBodyGetPosition(body1);
        pos2 = dBodyGetPosition(body2);
        //相対位置ベクトルを計算 box1->box2に向かう向きが正
        relPos[0] = pos2[0] - pos1[0];
        relPos[1] = pos2[1] - pos1[1];
        relPos[2] = pos2[2] - pos1[2];

        force1 = dBodyGetForce(body1);
        force2 = dBodyGetForce(body2);
        R1 = dBodyGetRotation(body1);
        R2 = dBodyGetRotation(body2);
        double theta1 = atan2(R1[4], R1[0]);
        double theta2 = atan2(R2[4], R2[0]);
        double delta_theta = theta2 -theta1;
            
        // printf("Xdistance=%f\n",abs(relPos[0]));
        // printf("Ydistance=%f\n",relPos[1]);
        double x1 = pos1[0] + 0.25;
        double x2 = pos2[0] - 0.25 * (cos(delta_theta) + sin(delta_theta));
        double y2 = sqrt(2) * 0.25 * sin(M_PI/4 - delta_theta);
        if (delta_theta <= M_PI/36 && -M_PI/36 <= delta_theta){
        //x方向の条件
            if (abs(relPos[0]) <= 1.0 && abs(relPos[1]) <= 1.0){
                if (force2[0] <= fx_lim && 0 <= force2[0]){
                dBodyAddForce(body1,force2[0]/2,0,0);
                dBodyAddForceAtPos(body2,-force2[0],force2[0]*0.25*sin(delta_theta),0, pos2[0] -0.25 * cos(theta2), pos2[1] -0.25 * sin(theta2),0.25);
                // fx2 = force2[0];
                }
                else if (force2[0] > fx_lim){
                dBodyAddForce(body1,force2[0]-fx_lim,0,0);
                    if (delta_theta >=0){
                        dBodyAddForceAtPos(body2,-fx_lim,fx_lim,0,pos2[0] -0.25 * cos(theta2), pos2[1] -0.25 * sin(theta2),0.25);
                    }
                    else if (delta_theta < 0){
                        dBodyAddForceAtPos(body2,-fx_lim,-fx_lim,0,pos2[0] -0.25 * cos(theta2), pos2[1] -0.25 * sin(theta2),0.25);
                    }
                }

            }
        //y方向の条件
            if (abs(relPos[0]) <= 1.0 && abs(relPos[1]) <= 1.0){
                if (abs(force2[1]) > 0){
                    if (force2[1] <= fy_lim &&  -fy_lim <= force2[1]){
                    // dBodyAddForceAtPos(body1,0,force2[1],0,pos1[0], pos1[1],0.5);
                    dBodyAddForce(body1,0,force2[1],0);
                    // fy2 = force2[1];
                    }
                    else if (force2[1] > fy_lim ){
                    dBodyAddForceAtPos(body1,0,force2[1] - fy_lim,0, pos1[0] +0.25 * cos(theta1), pos1[1] -0.25 * sin(theta1),0.25);
                    dBodyAddForceAtPos(body2,0,-fy_lim,0,pos2[0] -0.25 * cos(theta2), pos2[1] -0.25 * sin(theta2),0.25);
                    // fy2 = fy_lim;
                    }
                    else if (force2[1] < -fy_lim ){
                    dBodyAddForceAtPos(body1,0,fy_lim + force2[1],0, pos1[0] +0.25 * cos(theta1), pos1[1] -0.25 * sin(theta1),0.25);
                    dBodyAddForceAtPos(body2,0,fy_lim,0,pos2[0] -0.25 * cos(theta2), pos2[1] -0.25 * sin(theta2),0.25);
                    // fy2 = -fy_lim;
                    }
                }
            }
        }
        

        else if(delta_theta < -M_PI/36 && M_PI/36 < delta_theta){
        //このとき弱い粘着が発生する
            if (x2 - x1 < 1.0 && abs(relPos[1]) < 1.0){ 
                if (force2[0] <= fx_lim && 0 <= force2[0]){
                    if (delta_theta >= 0){
                        dBodyAddForceAtPos(body1,force2[0]/2,0,0,x1, pos2[1] + y2, 0.25);
                        dBodyAddForceAtPos(body2,-force2[0]/2, 0, 0, x2, pos2[1] + y2, 0.25);
                    }
                    else if (delta_theta < 0){
                        dBodyAddForceAtPos(body1,force2[0]/2,0,0,x1, pos2[1] - y2, 0.25);
                        dBodyAddForceAtPos(body2,-force2[0]/2, 0, 0, x2, pos2[1] - y2, 0.25);
                    }
                
                }
                else if (force2[0] > fx_lim){
                    if (delta_theta >= 0){
                        dBodyAddForceAtPos(body1,force2[0]-fx_lim,0,0,x1, pos1[1] , 0.25);
                        dBodyAddForceAtPos(body2,-fx_lim*(1.5 - relPos[0]),0,0, x2, pos2[1] + y2, 0.25);
                    }
                    if (delta_theta < 0){
                        dBodyAddForceAtPos(body1,force2[0]-fx_lim,0,0,x1, pos2[1] - y2, 0.25);
                        dBodyAddForceAtPos(body2,-fx_lim*(1.5 - relPos[0]),0,0, x2, pos2[1] - y2, 0.25);
                    }
                }
                if (abs(force2[1]) > 0){
                    if (force2[1] <= fy_lim &&  -fy_lim <= force2[1]){
                        if (delta_theta >= 0){
                        dBodyAddForceAtPos(body1,0,force2[1],0,x1, pos2[1] +y2, 0.25);
                        dBodyAddForceAtPos(body2,0,-force2[1],0,x2, pos2[1] +y2, 0.25);
                        }
                        else if (delta_theta < 0){
                        dBodyAddForceAtPos(body1,0,force2[1],0,x1, pos2[1] -y2, 0.25);
                        dBodyAddForceAtPos(body2,0,-force2[1],0,x2, pos2[1] -y2, 0.25);
                        }
                    }
                    else if (force2[1] > fy_lim ){
                        if (delta_theta >= 0){
                        dBodyAddForceAtPos(body1, 0, (force2[1] - fy_lim), 0, x1, pos2[1] + y2, 0.25); //なにかあったらy方向の力に（x2-x1）を掛ける
                        dBodyAddForceAtPos(body2,0,-fy_lim, 0, x2, pos2[1] + y2, 0.25);
                        }
                        else if (delta_theta < 0){
                        dBodyAddForceAtPos(body1, 0, (force2[1] - fy_lim), 0, x1, pos2[1] - y2, 0.25);
                        dBodyAddForceAtPos(body2,0,-fy_lim, 0, x2, pos2[1] - y2, 0.25);
                        }
                        
                    }
                    else if (force2[1] < -fy_lim ){
                        if (delta_theta >= 0){
                        dBodyAddForceAtPos(body1,0,fy_lim + force2[1],0, x1, pos2[1] + y2, 0.25);
                        dBodyAddForceAtPos(body2,0,fy_lim,0, x2, pos2[1] + y2, 0.25);
                        }
                        else if (delta_theta < 0){
                        dBodyAddForceAtPos(body1,0,fy_lim + force2[1],0, x1, pos2[1] - y2, 0.25);
                        dBodyAddForceAtPos(body2,0,fy_lim,0, x2, pos2[1] - y2, 0.25);
                        }
                    }
                }
            }
        }
    }

            

    

    double step(int action,dBodyID body1, dBodyID body2){
        // double x1 = state[0];
        // double x1_dot = state[1];
        // double y1 = state[2];
        // double y1_dot = state[3];
        // double theta1 = state[4];
        // double theta1_dot = state[5];
        // double x2 = state[6];
        // double x2_dot = state[7];
        // double y2 = state[8];
        // double y2_dot = state[9];
        // double theta2 = state[10];
        // double theta2_dot = state[11];
        double x1 = state[0];
        double x1_dot = state[1];
        double y1 = state[2];
        double y1_dot = state[3];
        double x2 = state[4];
        double x2_dot = state[5];
        double y2 = state[6];
        double y2_dot = state[7];
        double delta_theta = state[8];
        

        double dx = 2.0;
        double dy = 0.4;

        double forcex1 = 0;
        double forcex2 = 0;
        // forceL = dBodyGetForce(body1);
        // forceR = dBodyGetForce(body2);

        double fx_lim = 15.0;
        double fy_lim = 10.0;
        const dReal *bX1 = dBodyGetPosition(body1);
        const dReal *bX2 = dBodyGetPosition(body2);
        double fx = 0.0;
        double fy = 0.0;
        double posx, posy;

        if (action == 0){
            posx = 0, posy =0;
        }
        else if (action == 1){
            posx = bX2[0], posy = bX2[1] +0.25, fx = 10;
        }
        else if (action == 2){
            posx = bX2[0], posy = bX2[1] +0.25, fx = 10, fy =-10;
        }
        else if (action == 3){
            posx = bX2[0], posy = bX2[1] +0.25, fy = -10;
        }
        else if (action == 4){
            posx = bX2[0], posy = bX2[1] +0.25,fx = -10, fy = -10;
        }
        else if (action == 5){
            posx = bX2[0], posy = bX2[1] +0.25, fx = -10;
        }
        else if (action == 6){
            posx = bX2[0], posy = bX2[1] -0.25, fx = 10;
        }
        else if (action == 7){
            posx = bX2[0], posy = bX2[1] -0.25, fx = 10, fy = 10;
        }
        else if (action == 8){
            posx = bX2[0], posy = bX2[1] -0.25, fy = 10;
        }
        else if (action == 9){
            posx = bX2[0], posy = bX2[1] -0.25, fx = -10, fy = 10;
        }
        else if (action == 10){
            posx = bX2[0], posy = bX2[1] -0.25, fx = -10;
        }
        else if (action == 11){
            posx = bX2[0]+0.25, posy = bX2[1], fy = 10;
        }
        else if (action == 12){
            posx = bX2[0]+0.25, posy = bX2[1], fx = -10, fy = 10;
        }
        else if (action == 13){
            posx = bX2[0]+0.25, posy = bX2[1], fx = -10;
        }
        else if (action == 14){
            posx = bX2[0]+0.25, posy = bX2[1], fx = -10, fy = -10;
        }
        else if (action == 15){
            posx = bX2[0]+0.25, posy = bX2[1], fy = -10;
        }
        dBodyAddForceAtPos(body2, fx, fy, 0, posx, posy, 0.25);

        dSpaceCollide(space,0,&nearCallback);
        calculatestickness(body1, body2);
        dWorldStep(world,stepsize);
        dJointGroupEmpty(contactgroup);
        const dReal *X1 = dBodyGetPosition(body1);
        const dReal *X2 = dBodyGetPosition(body2);
        const dReal *V1 = dBodyGetLinearVel(body1);
        const dReal *V2 = dBodyGetLinearVel(body2);
        const dReal *R1  = dBodyGetRotation(body1);
        const dReal *R2  = dBodyGetRotation(body2);
        // const dReal *th1 = std::atan2(R[9], R[10]) ;
        // const dReal *the2;
        const dReal *omega1 = dBodyGetAngularVel(body1) ;
        const dReal *omega2 = dBodyGetAngularVel(body2);;
        
        x1 = X1[0];
        x1_dot = V1[0];
        // theta1 = std::atan2(R1[9], R1[10]);
        // theta1_dot = omega1[0];
        y1 = X1[1];
        y1_dot = V1[1];
        x2 = X2[0];
        x2_dot = V2[0];
        
        // theta2_dot = omega2[0];
        y2 = X2[1];
        y2_dot = V2[1];        
        double theta1 = std::atan2(R1[4], R1[0]);
        double theta2 = std::atan2(R2[4], R2[0]);

        delta_theta = theta2 - theta1;
        
        // state = {x1, x1_dot, theta1, theta1_dot, x2, x2_dot, theta2, theta2_dot};
        state = {x1, x1_dot,y1, y1_dot, x2, x2_dot,y2, y2_dot,delta_theta};

        bool terminated = terminate();
        double reward;
        // if (distance > 1.0){
        //     reward += 10;
        // }
        if (!terminated){
            reward = -5.0;
        }

        if (terminated){
            if (steps_beyond_terminated == 0){
                std::cerr <<"You are calling 'step()' even though this environment has already returned terminated = True. You should always call 'reset()' once you receive 'terminated' = True -- any further steps are undefined behavior."<<std::endl;
            }
            steps_beyond_terminated += 1;
            reward = 0;
        }
        return reward;
    }

    std::vector<double>getcurrentState(){
        return state;
    }

    std::vector<double> reset(){
        // state[0] = (rand() / (RAND_MAX + 1.0)) * 0.1 -0.05 - 0.5;
        // state[1] = (rand() / (RAND_MAX + 1.0)) * 0.1 -0.05;
        // state[2] = (rand() / (RAND_MAX + 1.0)) * 0.1 -0.05;
        // state[3] = (rand() / (RAND_MAX + 1.0)) * 0.1 -0.05;
        // state[4] = (rand() / (RAND_MAX + 1.0)) * 0.1 -0.05;
        // state[5] = (rand() / (RAND_MAX + 1.0)) * 0.1 -0.05;
        // state[6] = (-1.0 +(rand() / (RAND_MAX + 1.0)) *2.0) * 0.25;
        // state[7] = (rand() / (RAND_MAX + 1.0)) * 0.1 -0.05;
        // state[8] = 0;
        state[0] =  - 0.5;
        state[1] = 0 ;
        state[2] = 0;
        state[3] = 0;
        state[4] = 0;
        state[5] = 0;
        state[6] = 0;
        state[7] = 0;
        state[8] = 0;
        // state[9] = (rand() / (RAND_MAX + 1.0)) * 0.1 -0.05;
        // state[10] = (rand() / (RAND_MAX + 1.0)) * 0.1 -0.05;
        // state[11] = (rand() / (RAND_MAX + 1.0)) * 0.1 -0.05;
        

        // state ={ -0.5,0,0,0,0,0,0,0};
        // dBodySetPosition(box[0].body,state[0],0,0);
        // dBodySetPosition(box[1].body,state[4],0,0);
        steps_beyond_terminated = 0.0;
        return state;
    }
    
};
static void makeBox(dReal x0, dReal y0, dReal vx0, dReal vy0)
{
  dReal  z0 = 0.3;
  dMass mass;
  dMatrix3 R;

  box[0].body = dBodyCreate(world);
  dMassSetZero(&mass);
  dMassSetBoxTotal(&mass,m,side[0],side[1],side[2]);
  dBodySetMass(box[0].body,&mass);
//   dRFromAxisAndAngle(R, 0,1,0,0);
  dBodySetPosition(box[0].body, x0, y0, z0);
//   dBodySetRotation(box[0].body,R);
dBodySetLinearVel(box[0].body, vx0, vy0, 0);
  box[0].lx      = lx;
  box[0].ly      = ly;
  box[0].lz      = lz;
  box[0].geom   = dCreateBox(space,box[0].lx, box[0].ly, box[0].lz); // 球ジオメトリの生成
  dGeomSetBody(box[0].geom,box[0].body);         // ボディとジオメトリの関連付け
}

// 箱の生成2
static void makeBox2(dReal x0, dReal y0, dReal vx0, dReal vy0, dReal theta)
{
  dReal z0 = 0.3;
  dMass mass;
  dMatrix3 R2;
  
  box[1].body = dBodyCreate(world);
  dMassSetZero(&mass);
  dMassSetBoxTotal(&mass,m,side[0],side[1],side[2]);
  dBodySetMass(box[1].body,&mass);
  dRFromAxisAndAngle(R2, 0,0,1,theta);
  dBodySetPosition(box[1].body, x0, y0, z0);
  dBodySetRotation(box[1].body, R2);
  dBodySetLinearVel(box[1].body, vx0, vy0, 0);
  box[1].lx      = lx;
  box[1].ly      = ly;
  box[1].lz      = lz;
  box[1].geom   = dCreateBox(space,box[1].lx, box[1].ly, box[1].lz); // 球ジオメトリの生成
  dGeomSetBody(box[1].geom,box[1].body);         // ボディとジオメトリの関連付け
}
const double alpha = 0.2;
const double ganma = 0.99;

// Qテーブルを表すクラス
class QTable {
public:
    std::vector<std::vector<double>> table; //tableという２次元ベクトルを作る

    QTable(int num_states, int num_actions) : table(num_states, std::vector<double>(num_actions, 0.0)) {
        // std::mt19937 rng(std::time(0));
        // std::uniform_real_distribution<double> dist(-1.0, 1.0); // 例: -1.0から1.0の範囲でランダムな値を生成

        // // 各エントリをランダムな値で初期化
        // for (int i = 0; i < num_states; ++i) {
        //     for (int j = 0; j < num_actions; ++j) {
        //         table[i][j] = dist(rng);
        //     }
        // }
    
    }

    int argmax(const std::vector<double>& values) {
        return std::distance(values.begin(), std::max_element(values.begin(), values.end())); //argmaxを自作
    }

    int getAction(const std::vector<double>& state, int episode) {
        int nextAction;
        double epsilon = 0.5 * std::pow(0.99, episode);
        if (epsilon <= (rand() / static_cast<double>(RAND_MAX))) {
            nextAction = argmax(table[digitizeState(state)]);
        } else {
            nextAction = rand() % table[0].size(); //table[0].size は1行目の列数　すなわち２個
            std::cout<<"nextaction ="<<nextAction<<std::endl;
        }
        return nextAction;
    }

    void updateQValue(int state, int action, int nextState, int nextAction, double reward, double alpha, double ganma) { //Q値の更新はアクションの決定と分けた
        table[state][action] = (1 - alpha) * table[state][action] +
                               alpha * (reward + ganma * table[nextState][nextAction]);
    }

    int digitizeState(const std::vector<double>& state) {
        auto bins = [](double value, double clip_min, double clip_max, int num_bins) {
        if (value < clip_min) {
            return 0;
        } else if (value >= clip_max) {
            return num_bins - 1;
        } else {
            double bin_width = (clip_max - clip_min) / num_bins;
            return static_cast<int>((value - clip_min) / bin_width);
        }
        };

        int num_bins = 4;

        int digitized[9];
        // digitized[0] = bins(state[0], -2.4, 2.4, num_bins);
        // digitized[1] = bins(state[1], -3.0, 3.0, 2);
        // digitized[2] = bins(state[2], -2.4, 2.4, num_bins);
        // digitized[3] = bins(state[3], -3.0, 3.0, 2);
        // digitized[4] = bins(state[4],  0 , 90, 2);
        // digitized[5] = bins(state[5], -1.0, 1.0, 2);
        // digitized[6] = bins(state[6], -2.4, 2.4, num_bins);
        // digitized[7] = bins(state[7], -3.0, 3.0, 2);
        // digitized[8] = bins(state[8], -2.4, 2.4, num_bins);
        // digitized[9] = bins(state[9], -3.0, 3.0, 2);
        // digitized[10] = bins(state[10], 0, 90, 2);
        // digitized[11] = bins(state[11], -1.0, 1.0, 2);
        digitized[0] = bins(state[0], -3.0, 3.0, num_bins);
        // digitized[1] = bins(state[1], -3.0, 3.0, num_bins);
        digitized[1] = bins(state[1], -1.0, 3.0, 2);
        digitized[2] = bins(state[2], -3.0, 3.0, num_bins);
        // digitized[3] = bins(state[3], -3.0, 3.0, num_bins);
        digitized[3] = bins(state[3], -3.0, 3.0, 2);
        digitized[4] = bins(state[4],  -3.0, 3.0, num_bins);
        // digitized[5] = bins(state[5], -3.0, 3.0, num_bins);
        digitized[5] = bins(state[5], -3.0, 3.0, 2);
        digitized[6] = bins(state[6], -3.0, 3.0, num_bins);
        // digitized[7] = bins(state[7], -3.0, 3.0, num_bins);
        digitized[7] = bins(state[7], -3.0, 3.0, 2);
        digitized[8] = bins(state[8], -M_PI/12, M_PI/12, 3); //-15~+15度で3分割


        // 0~255に変換
        int result = 0;
        // for (int i = 0; i < 8; ++i) {
        //     result += digitized[i] * static_cast<int>(std::pow(4, i));
        // }
        result =  digitized[1]*4*4*4*4 + digitized[3] *4*4*4*4*2 + digitized[5]*4*4*4*4*2*2  + digitized[7]*4*4*4*4*2*2*2  + digitized[2]*4 + digitized[4]*4*4 + digitized[6] *4*4*4 +digitized[8]*4*4*4*4*2*2*2*2;
        return result;
    }
};

double calculateMean(const std::vector<double>& data){
    double sum = 0.0;
    for (const double& value : data){
        sum += value;
    }

    double mean = 0.0;
    if (!data.empty()){
        mean = sum / data.size();
    }
    return mean;
}
static void destroyBox0()
{
  dBodyDestroy(box[0].body);
  dGeomDestroy(box[0].geom);
}

//箱1の破壊
static void destroyBox1()
{
  dBodyDestroy(box[1].body);
  dGeomDestroy(box[1].geom);
}

int main() {
    const int num_states = 4 * 4 * 4 * 4 *  2* 2 * 2 * 2 *3;  // 仮の状態数
    const int num_actions = 16;  // アクション数    

  
    fp2 = fopen("Qtablexyt.csv","w");


    //ODEの設定
    dInitODE();
    world = dWorldCreate();
    dWorldSetGravity(world, 0,0,-9.8);
    space = dHashSpaceCreate(0);
    // contactgroup = dJointGroupCreate(0);
    ground = dCreatePlane(space, 0, 0, 1, 0);
    // makeBox();
    // makeBox2();
    

    // Qテーブルの初期化
    Stickness env;
    QTable qTable(num_states, num_actions);

    // エピソード関連のパラメータ
    const int num_episodes = 1000;
    const int max_number_of_steps = 2000;
    const int goal_average_steps = 200;
    const int num_consecutive_iterations = 100;
    const double alpha = 0.2;
    const double ganma = 0.99;
    std::vector<double> last_time_steps(num_consecutive_iterations, 0);
    std::vector<double> flags(num_consecutive_iterations, 0);
    std::vector<int> step_list;
    int episode =0;
    int success_flag =0;
    int final_flag = 0;
    
    double sumReward = 0.0;
    // エピソードのループ
    // for (int episode = 0; episode < num_episodes; ++episode) {
    while (true){       
        episode +=1; 
        // success_flag = 0;
        contactgroup = dJointGroupCreate(0);
        fp = fopen("xyt_episode_vs_step.csv","a");
        // 状態の取得
        std::vector<double> observation = env.reset();

        makeBox(observation[0],observation[2], observation[1],observation[4]);
        makeBox2(observation[4], observation[6],observation[5],observation[7],observation[8]);

        int state = qTable.digitizeState(observation);

        // アクションの取得
        int action = qTable.getAction(observation, episode);

        int episode_reward = 0;

        double prevXDistance = abs(observation[4] - observation[0]);
        double prevYDistance = abs(observation[6] - observation[2]);

        // タイムステップのループ
        for (int t = 0; t < max_number_of_steps; ++t) {
            // 報酬の取得
            double reward = env.step(action,box[0].body,box[1].body);
            observation = env.getcurrentState();
            double currentXDistance = abs(observation[4] - observation[0]);
            double currentYDistance = abs(observation[6] - observation[2]);

            if (env.terminate() ||t == max_number_of_steps - 1 ){
                printf("x1=%f,x2 = %f,xdistance = %f\n",observation[0],observation[4],abs(observation[4]-observation[0]));
                printf("y1=%f,y2 = %f,ydistance = %f\n",observation[2],observation[6],abs(observation[6]-observation[2]));
                if ((currentXDistance > 2.0)){
                    reward  +=100;
                    success_flag += 1;
                }
                else if (currentYDistance > 2.0){
                    reward += 50;
                    success_flag += 1;
                }
                else if ((currentXDistance <= 2.0) || (currentYDistance <= 2.0)){
                    reward = -500;
                    
                };
            }
            if (currentXDistance > prevXDistance){
                reward += 1.0;
            }
            if (currentYDistance > prevYDistance){
                reward += 1.0;
            }

                prevXDistance = currentXDistance;
                prevYDistance = currentYDistance;
            // 次のアクションの取得
            int nextstate = qTable.digitizeState(env.getcurrentState());
            int nextAction = qTable.getAction(env.getcurrentState(), episode);
            // std::cout << "next_action =" << nextAction << std::endl;
            // std::cout << "x2 =" << observation[6] << "x2 ="<< observation[0] << std::endl;
            

            // Qテーブルの更新
            qTable.updateQValue(state, action, nextstate, nextAction, reward, alpha, ganma);

            episode_reward += reward;
            //printf("state =%i, nextstate=%i, observation =%f %f %f %f %f %f %f %f, reward =%f,action = %i\n",state, nextstate,observation[0],observation[1],observation[2],observation[3],observation[4],observation[5],observation[6],observation[7],reward,nextAction);

            
            sumReward += std::pow(ganma, t+1)*reward;
            

            if (env.terminate() || t == max_number_of_steps - 1) {
                fp3 = fopen("xyt_sumreward.csv","a");
                std::cout << episode << " finished after " << t + 1 << " time steps / mean " << calculateMean(last_time_steps) << std::endl;
                std::cout << "Episode:" << episode << "  sum of reward =" << sumReward << "flags =" << calculateMean(flags)<<std::endl;
                 fprintf(fp3, "%i,%f\n",episode,sumReward);
                 sumReward = 0.0;
                // if (!sumRewards.empty()) {
                //     std::cout << "Episode:" << episode << ", sum of rewards:" << sumRewards.back() << std::endl;
                // } else {
                //     std::cout << "Episode:" << episode << ", sum of rewards: N/A" << std::endl;
                // }
                // fprintf(fp3, "%i,%f\n", episode, sumRewards.back());
                last_time_steps.erase(last_time_steps.begin());
                last_time_steps.push_back(t + 1);
                step_list.push_back(t + 1);
                flags.erase(flags.begin());
                flags.push_back(success_flag);
                
                

                printf("nextstate =%i, observation = \n",nextstate);
                dJointGroupDestroy(contactgroup);
                destroyBox0();
                destroyBox1();
                fprintf(fp, "%i,%i\n",episode,t+1);
                fclose(fp);
                fclose(fp3);
                // episode += 1;
                break;                
            }
            state = nextstate;
            action = nextAction;
            if (episode == 1){
                FILE *fp4 = fopen("actionlist_z.csv","a");
                fprintf(fp4, "%i\n",action);
                fclose(fp4);
            }
            if (episode == 100){
                FILE *fp5 = fopen("actionlist_h.csv","a");
                fprintf(fp5, "%i\n",action);
                fclose(fp5);
            }
            if (final_flag == 1){
                FILE *fp6 = fopen("actionlist_f.csv","a");
                fprintf(fp6, "%i\n",action);
                fclose(fp6);
                printf ("%i\n",action);
            }
            observation =env.getcurrentState();

        }
        // 成功判定（ここに成功判定の処理を追加）
        if (success_flag > 999){
            final_flag =1; 
        }
        if (success_flag >1000){
            std::cout << "Episode " << episode << " train agent successfully!" << std::endl;
            break;
        }
        // if (calculateMean(flags) >=0.99) {
        //     std::cout << "Episode " << episode << " train agent successfully!" << std::endl;
        //     break;
        // }

        // if (calculateMean(last_time_steps) >= goal_average_steps ) {
        //     std::cout << "Episode " << episode << " train agent successfully!" << std::endl;
        //     // for (int i = 0; i < num_states; i++){
        //     //     for (int j = 0; j < num_actions; j++){
        //     //         fprintf(fp2,"q[%i][%i] = %f\n", i,j,qTable.table[i][j]);
        //     //     }
        //     // }
        //     break;
        // }
        // if (episode == num_episodes -1){
        //     break;
        // }
    }
    for (int i = 0; i < num_states; i++){
                for (int j = 0; j < num_actions; j++){
                    fprintf(fp2,"%i,%i,%f\n", i,j,qTable.table[i][j]);
                }
            }
    dSpaceDestroy(space);
    dWorldDestroy(world);
    // fclose(fp);
    fclose(fp2);
    dCloseODE();
    return 0;
}
