//
//  main.cpp
//  algorithm
//
//  Created by Chen Xie on 2019/11/7.
//  Copyright © 2019 Chen Xie. All rights reserved.
//

#include <iostream>
#include <cmath>
int main(int argc, const char * argv[]) {
    // insert code here...
    double v[10];
    double r[10];
    double a[10];
    double M = 1;
    double w = 1;
    double H[10];
    double K[10];
    double V[10];
    double ast1 = 0, ast2 = 0;
    int i;
    double t = 0;
    double dt = 0.2;
    double step = 1000001;
    for(i = 0;i < 10;i++){
        v[i] = 5;
        r[i] = 0;
        a[i] = 0;
        H[i] = 0;
        K[i] = 0;
        V[i] = 0;
    }
    for (i = 0;i < step;i++){
        t = i * dt;
        // velocity verlet
        v[0] += 0.5 * dt * a[0];
        r[0] += dt * v[0];
        a[0] = -w * w * r[0];
        v[0] += 0.5 * dt * a[0];
        K[0] = 0.5 * M * v[0] *v[0];
        V[0] = 0.5 * M * w * w * r[0] * r[0];
        H[0] = K[0] + V[0];
        // original velocity verlet
        r[1] += dt * v[1] + 0.5 * dt * dt * a[1];
        v[1] += 0.5 * dt * a[1];
        a[1] = -w * w * r[1];
        v[1] += 0.5 * dt * a[1];
        K[1] = 0.5 * M * v[1] *v[1];
        V[1] = 0.5 * M * w * w * r[1] * r[1];
        H[1] = K[1] + V[1];
        if(i % 1 == 0){
     //   printf("%f, %f, %f,\n",t,H[0],H[1]);
        }
    }
    // leap frog initial condition
    double vst = 0;
    double rst3 = 0;
    a[2] = -w * w * r[2];
    v[2] -= 0.5 * a[2] * dt;
    // leap frog
    for(i =0; i < step;i++){
        t= i * dt;
        vst = v[2];
        rst3 = r[2];
        v[2] += a[2] * dt;
        r[2] += v[2] * dt;
        a[2] = -w * w * r[2];
        K[2] = 0.5 * M * (0.5 * vst + 0.5 * v[2]) * (0.5 * vst + 0.5 * v[2]);
        V[2] = 0.5 * M * w * w * rst3 * rst3;
        H[2] = K[2] + V[2];
        if(i % 1 == 0){
      //  printf("%f, %f,\n",t,H[2]);
        }
    }
    // Beeman initial condition
    a[3] = - w * w * r[3];
    ast1 = a[3];
    r[3] += -v[3] * dt + 0.5 * a[3] * dt * dt;
    a[3] = - w * w * r[3];
    ast2 = a[3];
    // Beeman
    for (i = 1; i < step;i++){
        t = i * dt;
        r[3] += v[3] * dt + (2.0/3.0 * ast1 - 1.0/6.0 * ast2) * dt * dt;
        a[3] = -w * w * r[3];
        v[3] += (1.0/3.0 * a[3] + 5.0/6.0 * ast1 - 1.0/6.0 * ast2) * dt;
        ast2 = ast1;
        ast1 = a[3];
        K[3] = 0.5 * M * v[3] * v[3];
        V[3] = 0.5 * M * w * w * r[3] * r[3];
        H[3] = K[3] + V[3];
        if(i % 1 == 0)  {
        printf("%f, %f,\n",t,H[3]);
        }
    }
    // Runge-Kutta-4
    double kr[4];
    double kv[4];
    for(i = 0; i< 4;i++){
        kr[i] = 0;
        kv[i] = 0;
    }
    for (i = 0;i < step;i++){
        t = i * dt;
        kr[0] = v[4];
        kv[0] = -r[4];
        kr[1] = v[4] + 0.5 * dt * kv[0];
        kv[1] = -(r[4] + 0.5 * dt * kr[0]);
        kr[2] = v[4] + 0.5 * dt * kv[1];
        kv[2] = -(r[4] + 0.5 * dt * kr[1]);
        kr[3] = v[4] + dt * kv[2];
        kv[3] = -(r[4] + dt * kr[2]);
        r[4] += dt * (kr[0] + 2 * kr[1] + 2 * kr[2] + kr[3]) / 6;
        v[4] += dt * (kv[0] + 2 * kv[1] + 2 * kv[2] + kv[3]) / 6;
        K[4] = 0.5 * M * v[4] * v[4];
        V[4] = 0.5 * M * w * w * r[4] * r[4];
        H[4] = K[4] + V[4];
        if(i % 100 == 0) {
    //    printf("%f, %f,\n",t,H[4]);
        }
    }
    // optimized leap frog initial condition
    double ast3 = 0;
    double ast4 = 0;
    double ast5 = 0;
    double vst2 = 0;
    double vst3 = 0;
    double rst4 = 0;
    double rst5 = 0;
    double rst6 = 0;
    rst4 = r[5];
    rst5 = r[5];
    rst4 -= r[5] - dt * v[5];
    rst5 += r[5] + dt * v[5];
    ast3 = -w * w * rst4;
    ast4 = -w * w * rst5;
    a[5] = -w * w * r[5];
    v[5] -= 0.5 * a[5] * dt;
    // optimized leap frog
    for (i = 0;i < step;i++){
        t = i * dt;
        vst2 = v[5];
        rst6 = r[5];
        ast5 = a[5];
        v[5] += dt * a[5] + (ast4 - 2 * a[5] + ast3) * dt /24;
        r[5] += dt * v[5] + (ast4 - a[5]) * dt * dt /24;
        a[5] =-w * w * r[5];
        K[5] = 0.5 * M * (0.5 * vst2 + 0.5 * v[5] - (ast4 - ast3) * dt/16) * (0.5 * vst2 + 0.5 * v[5] - (ast4 - ast3) * dt/16);
        V[5] = 0.5 * M * w * w * rst6 * rst6;
        H[5] = K[5] + V[5];
        ast3 = a[5];
        a[5] = ast4;
        vst3 = v[5] + 0.5 * ast5 * dt;
        rst5 = r[5] + vst3 * dt;
        ast4 = -w * w * rst5;
     //   printf("%f, %f\n",t,H[5]);
    }
    // Gear predictor-corrector
    double b1;
    double delta;
    double ap;
    b1 = -w * w * v[6];
    t = 0;
    for(i = 1;i < step;i++){
       t = i * dt;
       r[6] += v[6] * dt + 0.5 * a[6] * dt * dt + b1 * dt * dt * dt /6;
       v[6] += a[6] * dt + 0.5 * b1 * dt * dt;
       a[6] += b1 * dt;
       b1 = b1;
       ap = -w * w * r[6];
       delta = ap - a[6];
       r[6] += 19 * delta * dt * dt /240;
       v[6] += 3 * delta *dt /8;
       a[6] += delta;
       b1 += 1.5 * delta / dt;
       K[6] = 0.5 * M * v[6] * v[6];
       V[6] = 0.5 * M * w * w * r[6] * r[6];
       H[6] = K[6] + V[6];
        if(i % 100 == 0){
  //     printf("%f, %f,\n",t,H[6]);
        }
    }
    // 4th order Yoshida integrator
    return 0;
}
