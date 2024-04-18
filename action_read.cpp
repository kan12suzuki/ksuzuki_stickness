#include <ode/ode.h>
#include <drawstuff/drawstuff.h>
#include <iostream>
#include <cmath>
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#ifdef dDOUBLE                      // 単精度と倍精度の両方に対応する
#define dsDrawSphere dsDrawSphereD  // ためのおまじない
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
dReal lx = 0.5,ly = 0.5, lz = 0.5, m  = 1.0;
dsFunctions fn;
static dReal side[3] = {lx, ly, lz};
FILE *fp;
int t = 0;

typedef struct {       // MyObject構造体
  dBodyID body;        // ボディ(剛体)のID番号（動力学計算用）
  dGeomID geom;        // ジオメトリのID番号(衝突検出計算用）
  double lx, ly, lz, m;       // 縦，横，高さ，質量
} MyObject;
MyObject box[3];

typedef struct {        //箸のobject
  dBodyID body;
  dGeomID geom;
  double lxc, lyc, lzc, mc;
} Chopstick;
Chopstick chopstick;

// コールバック関数
static void nearCallback(void *data, dGeomID o1, dGeomID o2)
{
  static const int N = 10; // 接触点数の最大値
  dContact contact[N];     // 接触点

  int isGround = ((ground == o1) || (slope == o1))
               ||((ground == o2) || (slope == o2))
               ||((box[1].geom == o1) || (box[1].geom == o2))
               ||((box[0].geom == o1) || (box[0].geom == o2))
               ||((chopstick.geom == o1) || (chopstick.geom == o2));

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
  
}
//物体同士の距離を計算
dReal calculateDistance(dBodyID body1, dBodyID body2){
    const dReal*pos1 = dBodyGetPosition(body1);
    const dReal*pos2 = dBodyGetPosition(body2);
    
    dReal distance = sqrt(pow(pos2[0] - pos1[0], 2) + pow(pos2[1] - pos1[1], 2) + pow(pos2[2] - pos1[2], 2));
    return distance;
}
std::vector<int> readCSV(const std::string& filename) {
    std::ifstream file(filename); // ファイルを開く
    std::vector<int> data; // データを格納するための1次元配列

    if (file.is_open()) {
        std::string line;
        while (std::getline(file, line)) { // ファイルから1行ずつ読み取る
            std::istringstream iss(line);
            std::string token;
            while (std::getline(iss, token, ',')) { // カンマで区切られたトークンを取得
                int value = std::stoi(token); // 文字列を整数に変換
                data.push_back(value); // 整数を配列に追加
            }
        }
        file.close(); // ファイルを閉じる
    } else {
        std::cerr << "Failed to open file: " << filename << std::endl;
    }

    return data;
}
// 粘着力を計算
void calculatestickness(dBodyID body1, dBodyID body2){
  static dReal angle = 0;
  const dReal *pos,*R, *R1, *R2;
  dReal angle_deg;
  const dReal *force1,*force2,*force3,*T1,*T2,*T3,*Totalforce1,*Totalforce2;
  const dReal *pos1, *pos2,*pos3; 
  dReal relPos[3],relPos2[3];
  static dReal fx_val = 0, fy_val, fz_val ,k,fx_lim,fy_lim,fz_lim ,fx, fy, fz,fx2,fy2,fz2;

  // fx += -0.01;
  // if (fx_val<=5){fx_val += 0.01;}
  // else if (fx_val == 5){fx_val = 5;}
  // fx_val = 9.0;
  // fy_val = 0.0;
  // fz_val = 0.0;
  // k = 0.1;
  fx_lim = 20;
  fy_lim = 20;
  // printf("fx = %.3f\n",fx_val);
   
  pos1 = dBodyGetPosition(body1);
  pos2 = dBodyGetPosition(body2);
  //dBodyAddForceAtPos(box[1].body,25,0,0,pos2[0],pos2[1],0.5);
  //相対位置ベクトルを計算 box1->box2に向かう向きが正
  relPos[0] = pos2[0] - pos1[0];
  relPos[1] = pos2[1] - pos1[1];
  relPos[2] = pos2[2] - pos1[2];

  // dReal r = dSqrt(relPos[0]*relPos[0] + relPos[1]*relPos[1] + relPos[2]*relPos[2]);
  // fx = fx_val*relPos[0]/r;
  // fy = fy_val*relPos[1]/r;
  // fz = fz_val*relPos[2]/r;
  
  // dVector3 Impforce ;
  // dBodyAddForce(box[2].body,fx_val,fy_val,fz_val);
  //dWorldImpulseToForce(world,0.01,0.0,0.0,1.0, Impforce);
  //dBodyAddForceAtPos(box[1].body,Impforce[0],Impforce[1],Impforce[2],0,0.25,0);

  force1 = dBodyGetForce(body1);
  force2 = dBodyGetForce(body2);
  R1 = dBodyGetRotation(body1);
  R2 = dBodyGetRotation(body2);
  double theta1 = atan2(R1[4], R1[0]);
  double theta2 = atan2(R2[4], R2[0]);
  double delta_theta = theta2 -theta1;
  // fprintf(fp,"青にかかる力（前）=%f\n",force1[0]);
  // fprintf(fp,"緑にかかる力（前）=%f\n",force2[0]);
  
  dReal distance = calculateDistance(body1, body2); 
  // printf("Xdistance=%f\n",abs(relPos[0]));
  // printf("Ydistance=%f\n",relPos[1]);
  double x1 = pos1[0] + 0.25;
  double x2 = pos2[0] - 0.25 * (cos(delta_theta) + sin(delta_theta));
  double y2 = sqrt(2) * 0.25 * sin(M_PI/4 - delta_theta);

if (delta_theta <= M_PI/30 && -M_PI/30 <= delta_theta){
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
      // dBodyAddForce(body2,-fx_lim,0,0);
      //fx2 = fx_lim;
    }
    // if (force1[0] >= -fx_lim && force1[0] < 0){
    //    dBodyAddForce(body2,force1[0]/2,0,0);
    //    dBodyAddForce(body1,-force1[0]/2,0,0);     
    //   // fx = force1[0];
    // }
    // else if (force1[0] < -fx_lim){
    //   dBodyAddForce(body2, force1[0]+fx_lim,0,0);
    //   dBodyAddForce(body1, fx_lim*(1.5 - relPos[0]), 0,0);
    //   // fx = fx_lim;
    // }

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
    // else if (force2[1]<0 && -fy_lim < force2[1]){rr
    //   // dBodyAddRelForce(body1,0,force2[1],0);
    //   // dBodyAddRelForce(body2,0,-force2[1],0);
    //   fy2 = force2[1];
    // }
    else if (force2[1] < -fy_lim ){
      dBodyAddForceAtPos(body1,0,fy_lim + force2[1],0, pos1[0] +0.25 * cos(theta1), pos1[1] -0.25 * sin(theta1),0.25);
      dBodyAddForceAtPos(body2,0,fy_lim,0,pos2[0] -0.25 * cos(theta2), pos2[1] -0.25 * sin(theta2),0.25);
      // fy2 = -fy_lim;
    }
  }
  // else if (abs(force1[1]) > 0){
  //   if (force1[1] >= -fy_lim && force1[1] <= fy_lim){
  //     // dBodyAddForce(body1,0,-force1[1],0);
  //     // dBodyAddForce(body2,0,force1[1],0);
  //     dBodyAddForceAtPos(body2,0,force1[1]/2,0,pos2[0], pos[1],0.5);
  //     // fy = force1[1];
  //   }
  //   else if (force1[1] < -fy_lim){
  //     dBodyAddForce(body2, 0,(force1[1]+fy_lim),0);
  //     dBodyAddForce(body1, 0, fy_lim,0);
  //     // fy = -fy_lim;
  //   }
  //   else if (force1[1] > fy_lim ){
  //     dBodyAddForce(body2,0,(force1[1]-fy_lim),0);
  //     dBodyAddForce(body1,0,-fy_lim,0);
  //     // fy = fy_lim;
  //   }
  // }
  }
 }
  

 else if(delta_theta < -M_PI/30 && M_PI/30 < delta_theta){
  //このとき弱い粘着が発生する
 if (x2 - x1 < 1.0 && abs(relPos[1]) < 1.0){ 
    if (force2[0] <= fx_lim && 0 <= force2[0]){
      if (delta_theta >= 0){
        dBodyAddForceAtPos(body1,force2[0]/2,0,0,x1, pos2[1] + y2, 0.25);
        dBodyAddForceAtPos(body2,-force2[0]/2, 0, 0, x2, pos2[1] + y2, 0.25);
      }
      else if (delta_theta < 0){
        dBodyAddForceAtPos(body1,force2[0],0,0,x1, pos2[1] - y2, 0.25);
        dBodyAddForceAtPos(body2,-force2[0], 0, 0, x2, pos2[1] - y2, 0.25);
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
          dBodyAddForceAtPos(body1, 0, (force2[1] - fy_lim)/2, 0, x1, pos2[1] + y2, 0.25); //なにかあったらy方向の力に（x2-x1）を掛ける
          dBodyAddForceAtPos(body2,0,-fy_lim, 0, x2, pos2[1] + y2, 0.25);
        }
        else if (delta_theta < 0){
          dBodyAddForceAtPos(body1, 0, (force2[1] - fy_lim)/2, 0, x1, pos2[1] - y2, 0.25);
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
  
  
   Totalforce1 = dBodyGetForce(body1);
   Totalforce2 = dBodyGetForce(body2);
   T1 = dBodyGetTorque(body1);
   const dReal *X1,*X2,*V1,*V2;
   X1 = dBodyGetPosition(body1);
   X2 = dBodyGetPosition(body2);
   V1 = dBodyGetLinearVel(body1);
   V2 = dBodyGetLinearVel(body2);
  //  printf("物体にかかる力青: (%.3f, %.3f, %.3f)\n", force[0], force[1], force[2]);
  //  printf("物体にかかる力緑: (%.3f, %.3f, %.3f)\n", force2[0], force2[1], force2[2]);
  //  printf("物体にかかるトルク青: (%.3f, %.3f, %.3f)\n", T1[0], T1[1], T1[2]);
  // FILE *fp;
  // fp = fopen ("stickdate1.csv", "w");

  //  printf("物体にかかる力青: (%.3f, %.3f, %.3f)\n", Totalforce1[0], Totalforce1[1],Totalforce1[2]);
  //  printf("物体にかかる力緑: (%.3f, %.3f, %.3f)\n", Totalforce2[0], Totalforce2[1], Totalforce2[2]);
   printf( "物体の位置青:(%f,%f,%f)\n", X1[0], X1[1],X1[2]);
   printf( "物体の位置緑:(%f,%f,%f)\n", X2[0], X2[1],X2[2]);
  //  printf( "物体の速度青:(%f,%f,%f)\n", V1[0], V1[1],V1[2]);
  //  printf( "物体の速度緑:(%f,%f,%f)\n", V2[0], V2[1],V2[2]);
  //  std::cout<<"delta_theta = " << delta_theta * 180/M_PI <<std::endl;
   std::cout<<"dx = " << x2 -x1 <<std::endl;
  //  fprintf(fp,"物体にかかる力青: (%.3f, %.3f, %.3f)\n", Totalforce1[0], Totalforce1[1],Totalforce1[2]);
  //  fprintf(fp,"物体にかかる力緑: (%.3f, %.3f, %.3f)\n", Totalforce2[0], Totalforce2[1], Totalforce2[2]);
  

}
void ApplyForceToChopstick(dReal fx, dReal fy, dReal fz){
  dBodyAddForce(chopstick.body, fx, fy, fz);
}

void getaction(dBodyID body2){
    const dReal *bX2 = dBodyGetPosition(body2);
    double fx = 0.0;
    double fy = 0.0;
    double posx, posy;
    // std::vector<int> action_list = readCSV("actionlist_046z.csv");
    // std::vector<int> action_list = readCSV("actionlist_046h.csv");
    std::vector<int> action_list = readCSV("actionlist_f.csv");
      int action = action_list[t];
      printf("%i time steps, action = %i\n",t, action);

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
}

// シミュレーションループ
static void simLoop(int pause)
{
  if (!pause){
      dSpaceCollide(space,0,&nearCallback);  // 衝突検出関数
      getaction(box[1].body);
      calculatestickness(box[0].body,box[1].body);
      const dReal *X2 = dBodyGetPosition(box[1].body);
      dWorldStep(world,0.01);
      dJointGroupEmpty(contactgroup); // ジョイントグループを空にする
      t += 1;

  }
  dsSetColor(0.0,0.0,1.0);
  dsDrawBox(dBodyGetPosition(box[0].body),
                dBodyGetRotation(box[0].body), side);
 
  dsSetColor(0.0,1.0,0.0);
  dsDrawBox(dBodyGetPosition(box[1].body),
                dBodyGetRotation(box[1].body), side);
  // const dReal* rotation = dBodyGetRotation(box[1].body);
  //       double angleX = std::atan2(rotation[5], rotation[8]) * (180.0 / M_PI);
  //       double angleY = std::asin(-rotation[2]) * (180.0 / M_PI);
  //       double angleZ = std::atan2(rotation[1], rotation[0]) * (180.0 / M_PI);

  //       std::cout << "物体の角度(XYZ): " << angleX << ", " << angleY << ", " << angleZ << std::endl;
        
        dMatrix3 R;
        dRFromAxisAndAngle(R, 1,0,0,M_PI/3);
        // printf("%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",R[0],R[1],R[2],R[3],R[4],R[5],R[6],R[7],R[8],R[9],R[10],R[11]);
        double angleX = std::atan2(R[9], R[10]) * (180.0 / M_PI);
        // printf("%f",angleX);

//   dsSetColor(1.2, 1.2, 0.0);
//   dVector3 sides;
//   dGeomBoxGetLengths(slope, sides);
//   dsDrawBox(dGeomGetPosition(slope), dGeomGetRotation(slope), sides);

}

// 箱の生成
static void makeBox()
{
  dReal x0 =-0.5, y0 = 0, z0 = 0.3;
  dMass mass;
  dMatrix3 R2;

  box[0].body = dBodyCreate(world);
  dMassSetZero(&mass);
  dMassSetBoxTotal(&mass,m,side[0],side[1],side[2]);
  dBodySetMass(box[0].body,&mass);
//   dRFromAxisAndAngle(R2, 0,1,0,0);
  dBodySetPosition(box[0].body, x0, y0, z0);
//   dBodySetRotation(box[0].body,R2);
  box[0].lx      = lx;
  box[0].ly      = ly;
  box[0].lz      = lz;
  box[0].geom   = dCreateBox(space,box[0].lx, box[0].ly, box[0].lz); // 球ジオメトリの生成
  dGeomSetBody(box[0].geom,box[0].body);         // ボディとジオメトリの関連付け
}

// 箱の生成2
static void makeBox2()
{
  std::srand(static_cast<unsigned int>(std::time(nullptr)));
  dReal x0 = 0  , y0 = 0, z0 = 0.3;
  dMass mass;
  dMatrix3 R;
  dReal ang =  static_cast<double>(std::rand())/RAND_MAX * (M_PI/30);
  
  box[1].body = dBodyCreate(world);
  dMassSetZero(&mass);
  dMassSetBoxTotal(&mass,m,side[0],side[1],side[2]);
  dBodySetMass(box[1].body,&mass);
  dRFromAxisAndAngle(R, 0,0,1,ang);
  fprintf(fp,"%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",R[0],R[1],R[2],R[3],R[4],R[5],R[6],R[7],R[8],R[9],R[10],R[11]);
  dBodySetPosition(box[1].body, x0, y0, z0);
  dBodySetRotation(box[1].body, R);
  box[1].lx      = lx;
  box[1].ly      = ly;
  box[1].lz      = lz;
  box[1].geom   = dCreateBox(space,box[1].lx, box[1].ly, box[1].lz); // 球ジオメトリの生成
  dGeomSetBody(box[1].geom,box[1].body);         // ボディとジオメトリの関連付け
}
//箱0の破壊
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


//　斜面の生成
static void makeslope()
{
  dReal sx = 100, sy = 2, sz = 0.01;  //絶対座標系に沿った長さ
  dReal x = 0, y = 0, z = 0;          //重心の座標
  dReal ax = 0, ay = 1, az = 0;       // 回転軸ベクトル
  dReal angle = 45.0 * M_PI / 180.0;  // 回転角
  
  slope = dCreateBox(space, sx, sy, sz);
  dMatrix3 R;
  dRFromAxisAndAngle(R, ax, ay, az, angle);
  dGeomSetPosition(slope, x, y, z);
  dGeomSetRotation(slope, R);

}
void restart(){
  dJointGroupDestroy(contactgroup);
  destroyBox0();
  destroyBox1();

  contactgroup = dJointGroupCreate(0);
  makeBox();
  makeBox2();
}
//コマンド関数
void command(int cmd)
{

  float xyz[3], hpr[3];
  const dReal *X1 = dBodyGetPosition(box[0].body);
  const dReal *X2 = dBodyGetPosition(box[1].body);
  switch(cmd){
    // case 'd':dBodyAddForceAtPos(box[0].body,30,0,-20,X1[0],X1[1],0.5)  ;break;
    // case 'a':dBodyAddForceAtPos(box[0].body, -30,0,-20,X1[0],X1[1],0.5);break;
    // case 'w':dBodyAddForceAtPos(box[0].body, 0,20,-20,X1[0],X1[1],0.5) ;break;
    // case 's':dBodyAddForceAtPos(box[0].body, 0,-20,-20,X1[0],X1[1],0.5);break;
    // case 'q':dBodyAddForce(box[0].body, 0,0,20) ;break;
    // case 'e':dBodyAddForce(box[0].body, 0,0,-20);break;
    // // case 'l':dBodyAddForce(box[1].body,15,0,0)  ;break;
    // // case 'j':dBodyAddForce(box[1].body, -10,0,0);break;
    // // case 'i':dBodyAddForce(box[1].body, 0,10,0) ;break;
    // // case 'k':dBodyAddForce(box[1].body, 0,-10,0);break;
    case 'l':dBodyAddForceAtPos(box[1].body,30,0,0,X2[0],X2[1],0.25)  ;break;
    case 'j':dBodyAddForceAtPos(box[1].body, -10,0,0,X2[0],X2[1] -0.25,0.25);break;
    // case ';':dBodyAddForceAtPos(box[1].body,10,0,0,X2[0],X2[1] -0.25,0.25)  ;break;
    // case 'h':dBodyAddForceAtPos(box[1].body, -10,-10,0,X2[0],X2[1] +0.25,0.25);break;
    case 'i':dBodyAddForceAtPos(box[1].body, 0,10,0,X2[0],X2[1] -0.25,0.25) ;break;
    case 'k':dBodyAddForceAtPos(box[1].body, 0,-10,0,X2[0],X2[1] - 0.25,0.25);break;
    // case 'u':dBodyAddRelForce(box[1].body, 0,0,20);break;
    // case 'o':dBodyAddForceAtPos(box[1].body, 0,0,20,0.25,0.25,0.25);break;
    // case 't':dBodyAddForce(box[1].body, 10,10,0);break;
    case 'r':restart()                          ;break;
    case '1':dBodyAddForceAtPos(box[1].body,10,0,0,X2[0],X2[1] +0.25,0.25)  ;break;
    case '2':dBodyAddForceAtPos(box[1].body,10,-10,0,X2[0],X2[1] +0.25,0.25)  ;break;
    case '3':dBodyAddForceAtPos(box[1].body,0,-10,0,X2[0],X2[1] +0.25,0.25)  ;break;
    case '4':dBodyAddForceAtPos(box[1].body,-10,-10,0,X2[0],X2[1] +0.25,0.25)  ;break;
    case '5':dBodyAddForceAtPos(box[1].body,-10,0,0,X2[0],X2[1] +0.25,0.25)  ;break;
    case '6':dBodyAddForceAtPos(box[1].body,10,0,0,X2[0],X2[1] -0.25,0.25)  ;break;
    case '7':dBodyAddForceAtPos(box[1].body,10,10,0,X2[0],X2[1] -0.25,0.25)  ;break;
    case '8':dBodyAddForceAtPos(box[1].body,0,10,0,X2[0],X2[1] -0.25,0.25)  ;break;
    case '9':dBodyAddForceAtPos(box[1].body,-10,10,0,X2[0],X2[1] -0.25,0.25)  ;break;
    case '0':dBodyAddForceAtPos(box[1].body,-10,0,0,X2[0],X2[1] -0.25,0.25)  ;break;
    case 'q':dBodyAddForceAtPos(box[1].body,0,10,0,X2[0]+0.25,X2[1] ,0.25)  ;break;
    case 'w':dBodyAddForceAtPos(box[1].body,-10,10,0,X2[0]+0.25,X2[1] ,0.25)  ;break;
    case 'e':dBodyAddForceAtPos(box[1].body,-10,0,0,X2[0]+0.25,X2[1] ,0.25)  ;break;
    case 't':dBodyAddForceAtPos(box[1].body,-10,-10,0,X2[0]+0.25,X2[1] ,0.25)  ;break;
    case 'y':dBodyAddForceAtPos(box[1].body,0,-10,0,X2[0]+0.25,X2[1] ,0.25)  ;break;
    // case 'w':dBodyAddForceAtPos(box[1].body,-10,10,0,X2[0]+0.25,X2[1] ,0.25)  ;break;
    case 'v':dsGetViewpoint(xyz,hpr);
             printf("xyz=%4.2f %4.2f %4.2f", xyz[0],xyz[1],xyz[2]);
             printf("hpr=%6.2f %6.2f %6.2f", hpr[0],hpr[1],hpr[2]);break;
    // case 'd':dBodyAddForceAtPos(box[0].body,30,0,-20,X1[0],X1[1],0.5)  ;break;
                 
  }
}
void start()                                  /*** 前処理　***/
{
  // static float xyz[3] = {-8, 0, 10};
  // static float hpr[3] = {-180, -90,0};
  
  static float xyz[3] = {0.26,0.0,7};         // 視点の位置
  static float hpr[3] = {90, -90, 0};          // 視線の方向

  // static float xyz[3] = {-0.3,-3.66,0.5};         // 視点の位置
  // static float hpr[3] = {90, -0, 0};
  dsSetViewpoint(xyz,hpr);                     // カメラの設定
}



void setDrawStuff()           /*** 描画関数の設定 ***/
{
  fn.version = DS_VERSION;    // ドロースタッフのバージョン
  fn.start   = &start;        // 前処理 start関数のポインタ
  fn.step    = &simLoop;      // simLoop関数のポインタ
  fn.command = &command;
  fn.path_to_textures = "../../drawstuff/textures"; // テクスチャ
}

int main(int argc, char **argv)
{
  setDrawStuff();

  dInitODE(); // ODEの初期化
  world = dWorldCreate();
  dWorldSetGravity(world,0,0,-9.81);
  dMatrix3 R;
  
  fp = fopen ("stickdate1.csv", "w");
  dRFromAxisAndAngle(R, 0,0,1,M_PI/3);
  fprintf(fp,"%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",R[0],R[1],R[2],R[3],R[4],R[5],R[6],R[7],R[8],R[9],R[10],R[11]);

  space        = dHashSpaceCreate(0);   // 衝突用空間の創造
  contactgroup = dJointGroupCreate(0);  // ジョイントグループの生成
  ground = dCreatePlane(space,0,0,1,0); // 平面ジオメトリの生成
  //makeslope();
  makeBox(); 			// 球の作成
  makeBox2();

  dsSimulationLoop(argc,argv,920, 720,&fn);
  dSpaceDestroy(space);
  dWorldDestroy(world);
  fclose(fp);
  dCloseODE(); // ODEの終了
  return 0;
}
