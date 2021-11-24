#include "FreeSpace.hpp"

int filter_x[28]={-1,0,1,-3,-2,2,3,-4,4,-4,4,-5,5,-5,5,-5,5,-1,0,1,-3,-2,2,3,-4,4,-4,4};
int filter_y[28]={-5,-5,-5,-4,-4,-4,-4,-3,-3,-2,-2,-1,-1,0,0,1,1,5,5,5,4,4,4,4,3,3,2,2};
int all_x[89]={-1,0,1, 
                -3,-2,-1,0,1,2,3, 
                -4,-3,-2,-1,0,1,2,3,4, 
                -4,-3,-2,-1,0,1,2,3,4, 
                -5,-4,-3,-2,-1,0,1,2,3,4,5,
                -5,-4,-3,-2,-1,0,1,2,3,4,5,
                -5,-4,-3,-2,-1,0,1,2,3,4,5,
                -1,0,1,
                -3,-2-1,0,1,2,3,
                -4,-3,-2,-1,0,1,2,3,4,
                -4,-3,-2,-1,0,1,2,3,4};
int all_y[89]={-5,-5,-5,
                -4,-4,-4,-4,-4,-4,-4,
                -3,-3,-3,-3,-3,-3,-3,-3,-3,
                -2,-2,-2,-2,-2,-2,-2,-2,-2,
                -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
                0,0,0,0,0,0,0,0,0,0,0,
                1,1,1,1,1,1,1,1,1,1,1,
                5,5,5,
                4,4,4,4,4,4,4,
                3,3,3,3,3,3,3,3,3,
                2,2,2,2,2,2,2,2,2};

LivoxFreeSpace::LivoxFreeSpace()
{
    this->pVImg=(unsigned char*)calloc(DN_SAMPLE_IMG_NX*DN_SAMPLE_IMG_NY*DN_SAMPLE_IMG_NZ,sizeof(unsigned char));
}
LivoxFreeSpace::~LivoxFreeSpace()
{
    if(this->pVImg!=NULL)
    {
        free(this->pVImg);
    }
}


//                   输入             分辨率为1度，半径为50m的圆     一个float的vector用于存放最终结果
void LivoxFreeSpace::FreeSpaceFilter(float* free_space_small, int n , std::vector<float> & free_space)
{
    clock_t t0, t1, t2, t3, t4;
    t0 = clock();
    float pixel_size = 0.2, delta_d_in_r = 0.13, delta_r = 0.15; //delta_d_in_r is smaller than pixel_size and delta_r, to make sure all pixels are covered
    
    // 大小为 100 * 100m， 分辨率为0.2m的像素栅格，用来标记哪些是占据，哪些是free
    Eigen::MatrixXi src = Eigen::MatrixXi::Zero(100/pixel_size, 100/pixel_size); 
    Eigen::MatrixXi dst = Eigen::MatrixXi::Zero(100/pixel_size, 100/pixel_size);
    
    // 1度的弧长会随着半径的增大而增大，但是像素比例为固定大小0.2，因此半径越大时需要填充的像素就越多
    // 根据 l = alpha * r 可以知道，弧长为固定0.2时，对应的delta角度为 alpha = 0.2 / r
    // 也就是当半径为r时，扇形弧长变大，需要以一定的delta角度，逐个判断栅格是否占据，  **注意：这里的delta是弧度**
    std::vector<float> delta_t;
    for (float j = 0.0001; j < 50; j += delta_r) // Prepare the delta theta of different radius
    {
        delta_t.push_back(delta_d_in_r/j); // delta_d_in_r 小于 0.2 是为了保险起见，保证没有栅格漏掉
    } 
    for (int i = 0; i < 360; i++)
    {
        // 与左右两边相邻的扇形的r进行比较，选较小值 
        float r = min(free_space_small[i], free_space_small[(i + 1) % n]);
        r = min(r, free_space_small[(i - 1 + n) % n]);
        r = sqrt(r);


        int k = 0;
        // 在半径上按照 delta_r = 0.15 的步长去访问每一个栅格， free为1
        for (float j = 0; j < r - 0.5; j += delta_r)
        {
            float dt = delta_t[k++];
            float theta = (i - 180)*FREE_PI/180.0; // 角度修改为-pi~pi
            for (float t = theta - 0.01; t < theta + 0.01; t+=dt)
            {
                float x = j*cos(t);
                float y = j*sin(t);
                int m =int((50.0 - x) / pixel_size);
                int nn =int((50.0 - y) / pixel_size);
                src(m, nn) = 1;
            }
        }
    }

    t1 = clock();
    for (int i = 0; i < 360; i++)
    {
        for (float j = 0; j < 49; j += delta_r)  // 这两个for循环是为了遍历每一个角度的每一个半径的梯形部分
        {  
            // 计算角度为 i， 半径长度为 j 的栅格在 Matrix 中的坐标
            float x = j * cos((i - 180)*FREE_PI/180.0);
            float y = j * sin((i - 180)*FREE_PI/180.0);
            int m =int((50.0 - x) / pixel_size);
            int nn =int((50.0 - y) / pixel_size);

            // 计算该栅格所在的角度
            int theta = int(atan2f(y, x) * 180.0 / FREE_PI + 180.0 + 0.5); // 角度，四舍五入取整
            theta = theta % n; // 360的index为0

            // 找相邻扇形的较小r
            float r = min(free_space_small[theta], free_space_small[(theta + 1) % n]);
            r = min(r, free_space_small[(theta - 1 + n) % n]);

            if (r > j*j + 1)  // +1的原因是只对1m以外的栅格进行处理 
            {
                // 判断周边栅格是否为free，大约是一个1m的圆，28个大概相当于这个圆的轮廓
                int result = 0;
                for (int k = 0; k < 28; k++)  
                {
                    result += src(m + filter_x[k], nn + filter_y[k]);
                }
                if (result < 28) // check if position (m, nn) is in free space
                    break;
                
                // 在dst矩阵中对这个半径为1m的圆内的所有点进行判断，全部置为free状态， 相当于降噪功能
                for (int k = 0; k < 89; k++) 
                {
                    dst(m+all_x[k], nn+all_y[k]) = max(1, dst(m+all_x[k], nn+all_y[k]));
                }
                dst(m, nn) = 2;
            }
        }
    }


    t2 = clock();

    for (int i = 0; i < dst.rows(); i++)
    {
        for (int j = 0; j < dst.cols(); j++)
        {
            if (dst(i, j) > 0)
            {
                float x = (100.0 - i*pixel_size) - 50.0;
                float y = (100.0 - j*pixel_size) - 50.0;
                free_space.push_back(x);
                free_space.push_back(y);
                free_space.push_back(255);                    
            }
        }
    }
    t3 = clock();
    // printf("filter time: %f, generate map: %f, conv: %f, fs generate: %f\n\n", 1000.0*(t3 - t0) / CLOCKS_PER_SEC,
    //         1000.0*(t1 - t0) / CLOCKS_PER_SEC, 1000.0*(t2 - t1) / CLOCKS_PER_SEC, 1000.0*(t3 - t2) / CLOCKS_PER_SEC);
}


//                            输入   非地面点  非地面点数量  大小为360的float数组  float数组长度
void LivoxFreeSpace::FreeSpace(float* fPoints, int n, float* free_space, int free_space_n)
{
    int thetaId;
    float distance;

    // 输出结果初始化  一个分辨率为1度，半径为50m的圆
    for(int ii=0; ii < free_space_n; ii++)
    {
        free_space[ii] = 2500; 
    }

    // 遍历所有的非地面点
    for(int pid=0;pid<n;pid++)
    {
        // z值小于3m
        if(fPoints[pid*4+2] < 3) // points of high tree, buildings are rejected
        {
            // 主要是为了防止有一些点是打在车身上的
            if (abs(fPoints[pid*4 + 1]) < 1.2 && abs(fPoints[pid*4]) < 2.5) // reject the near points of robot
                continue;
            distance = fPoints[pid*4]*fPoints[pid*4] + fPoints[pid*4+1]*fPoints[pid*4+1];

            // 角度范围是逆时针从0到2pi  然后转换为角度，  +0.5的操作是因为int是向下取整，去掉尾巴的方式，而此处希望四舍五入取整
            thetaId = int((atan2f(fPoints[pid*4+1], fPoints[pid*4]) + FREE_PI) * 180.0 / FREE_PI + 0.5);
            thetaId = thetaId % free_space_n; // 360的index为0

            // 找该区间（角度一度）内距离最近的一个点， 如果该区间没有点则最远距离为50m
            if(free_space[thetaId] > distance && distance > 1)
            {
                free_space[thetaId] = distance;
            }
        }
    }

}

int LivoxFreeSpace::GenerateFreeSpace(float* fPoints1, int pointNum, std::vector<float> & free_space)
{
    clock_t t0, t1, t2, t3, t4;
    t0 = clock();
    // down sampling
    float *fPoints2=(float*)calloc(pointNum*4,sizeof(float));

    int *idtrans1=(int*)calloc(pointNum,sizeof(int));
    int *idtransx=(int*)calloc(pointNum,sizeof(int));
    int *idtransy=(int*)calloc(pointNum,sizeof(int));
    int *idtrans2=(int*)calloc(pointNum,sizeof(int));

    int pntNum = 0;

    //                                 500 * 200 * 100
    this->pVImg=(unsigned char*)calloc(DN_SAMPLE_IMG_NX*DN_SAMPLE_IMG_NY*DN_SAMPLE_IMG_NZ,sizeof(unsigned char));
    //                     200 * 80 * 100
    std::vector<int> count(DENOISE_IMG_NX*DENOISE_IMG_NY*DENOISE_IMG_NZ, 0);

    // 遍历输入点云， 划分颗粒度较大的栅格，为后面的降噪做准备
    for(int pid=0;pid<pointNum;pid++)
    {
        // 将每一个点划分到 -50 < x < 150m 分辨率1m    -40 < y < 40m 分辨率1m   -1 < z < 10m 分辨率0.2m的栅格中
        int ix=(fPoints1[pid*4]+DN_SAMPLE_IMG_OFFX)/DENOISE_IMG_DX; 
        int iy=(fPoints1[pid*4+1]+DN_SAMPLE_IMG_OFFY)/DENOISE_IMG_DY; 
        int iz=(fPoints1[pid*4+2]+DN_SAMPLE_IMG_OFFZ)/DENOISE_IMG_DZ;
        idtrans1[pid]=-1;
        if((ix>=0)&&(ix<DENOISE_IMG_NX)&&(iy>=0)&&(iy<DENOISE_IMG_NY)&&(iz>=0)&&(iz<DENOISE_IMG_NZ)) 
        {
            // 计算其在三维栅格中的索引
            int idx = iz*DENOISE_IMG_NX*DENOISE_IMG_NY+iy*DENOISE_IMG_NX+ix;

            // 根据每一个点的原始索引 记录其在三维栅格中的索引
            idtrans1[pid]=idx;

            // 统计每一个三维栅格中点的个数
            count[idx]++;     
        }
    }

    // 遍历输入点云, 降噪， 去掉稀疏点
    for(int pid=0;pid<pointNum;pid++)
    {
        // 在输入的点云中， 若其属于划定的三维栅格 且 该栅格内的点的个数小于3个   则将其重置到原点处
        if(idtrans1[pid] > -1 && count[idtrans1[pid]] < 3)
        {
            fPoints1[pid*4] = 0;
            fPoints1[pid*4 + 1] = 0;
            fPoints1[pid*4 + 2] = 0;

        }
    }

    // 遍历输入点云
    for(int pid=0;pid<pointNum;pid++)
    {
        // 将每一个点划分到 -50 < x < 150m 分辨率0.4m    -40 < y < 40m 分辨率0.4m   -1 < z < 10m 分辨率0.2m的栅格中
        int ix=(fPoints1[pid*4]+DN_SAMPLE_IMG_OFFX)/DN_SAMPLE_IMG_DX;
        int iy=(fPoints1[pid*4+1]+DN_SAMPLE_IMG_OFFY)/DN_SAMPLE_IMG_DY;
        int iz=(fPoints1[pid*4+2]+DN_SAMPLE_IMG_OFFZ)/DN_SAMPLE_IMG_DZ;

        idtrans1[pid]=-1;
        if((ix>=0)&&(ix<DN_SAMPLE_IMG_NX)&&(iy>=0)&&(iy<DN_SAMPLE_IMG_NY)&&(iz>=0)&&(iz<DN_SAMPLE_IMG_NZ))
        {
            // 存放在栅格图中的索引， 此时已变成细粒度的栅格图
            idtrans1[pid] = iz*DN_SAMPLE_IMG_NX*DN_SAMPLE_IMG_NY+iy*DN_SAMPLE_IMG_NX+ix;

            // 根据每一个点的原始索引 记录其在更细粒度的三维栅格中x方向的索引
            idtransx[pid] = ix;

            // 根据每一个点的原始索引 记录其在更细粒度的三维栅格中y方向的索引
            idtransy[pid] = iy;

            // pVImg是生成的细粒度三维栅格图    fPoints2只存落在栅格中的第一个点，其数量等于栅格中有点的栅格数，
            if(pVImg[idtrans1[pid]] == 0)//没有访问过，肯定栅格内会有重复的，所以fPoints2只取第一个
            {
                fPoints2[pntNum*4] = fPoints1[pid*4];
                fPoints2[pntNum*4 + 1] = fPoints1[pid*4+1];
                fPoints2[pntNum*4 + 2] = fPoints1[pid*4+2];
                fPoints2[pntNum*4 + 3] = fPoints1[pid*4+3];

                // pntNum初始时为0，在栅格第一次被访问时会加1
                // idtrans2 存 fPoints2 中对应位置的点在输入点云中的索引
                idtrans2[pntNum] = pid;
                pntNum++; 
            }
            
            // 最终只要栅格中有点，则在 pVImg 的对应位置都会被置1
            pVImg[idtrans1[pid]] = 1;
        }
    }

    t1 = clock();

    // 进行地面点分割             pntNum 表示有点的栅格的数量
    int *pLabelGnd=(int*)calloc(pntNum,sizeof(int));
    int ground_point_num = GroundSegment(pLabelGnd, fPoints2, pntNum, 1.0);

    t2 = clock();

    // agnum 表示非地面点的数量    fPoints3 存非地面点
    int agnum = pntNum - ground_point_num;
    float *fPoints3 = (float*)calloc(agnum*4,sizeof(float));
    int agcnt=0;
    for(int ii=0;ii<pntNum;ii++)
    {
        if(pLabelGnd[ii]==0)
        {
            fPoints3[agcnt*4]=fPoints2[ii*4];
            fPoints3[agcnt*4+1]=fPoints2[ii*4+1];
            fPoints3[agcnt*4+2]=fPoints2[ii*4+2];
            fPoints3[agcnt*4+3]=fPoints2[ii*4+3];
            agcnt++;
        }
        
    }
    float *free_space_small = (float*)calloc(360,sizeof(float));
    // 输入          非地面点  非地面点数量  大小为360的float数组  float数组长度
    this->FreeSpace(fPoints3, agnum, free_space_small, 360);

    //                   前一步得到的分辨率为1度，半径为50m的圆     一个float的vector用于存放最终结果
    this->FreeSpaceFilter(free_space_small, 360, free_space);

    free(fPoints2);
    free(idtrans1);
    free(idtrans2);
    free(idtransx);
    free(idtransy);
    free(fPoints3);
    free(pLabelGnd);
    free(this->pVImg);
    free(free_space_small);
    std::vector<int>().swap(count);
    t3 = clock();
    // printf("FreeSpace total Time: %f, Downsample: %f, Ground Segment: %f, FreeSpace: %f\n\n", 1000.0*(t3 - t0) / CLOCKS_PER_SEC, 
    //         1000.0*(t1 - t0) / CLOCKS_PER_SEC, 1000.0*(t2 - t1) / CLOCKS_PER_SEC, 1000.0*(t3 - t2) / CLOCKS_PER_SEC);
}

/*
int LivoxFreeSpace::GroundSegment(int* pLabel,float *fPoints,int pointNum,float fSearchRadius)
Fast ground segmentation using rule-based & plane fitting method 
*/

//                                pLabel 要输出的结果     fPoints 有点的栅格中第一个点   pointNum 有点的栅格的数量
int LivoxFreeSpace::GroundSegment(int* pLabel,float *fPoints,int pointNum,float fSearchRadius)
{
    int gnum=0;

    //                                 24 * 20
    float *pGndImg1 = (float*)calloc(GND_IMG_NX1*GND_IMG_NY1,sizeof(float));
    int *tmpLabel1 = (int*)calloc(pointNum,sizeof(int));
    
    // 初始化 pGndImg1 默认全部置为100
    for(int ii = 0; ii < GND_IMG_NX1 * GND_IMG_NY1; ii++)
    {
        pGndImg1[ii]=100;
    }

    // 遍历输入进来的每一个点 
    for(int pid = 0;pid < pointNum; pid++)
    {
        // 将每一个点划分到   -40 < x < 56m 分辨率4m     -40 < y < 40m 分辨率4m   的2D栅格中
        int ix= (fPoints[pid*4]+GND_IMG_OFFX1)/(GND_IMG_DX1+0.000001);
        int iy= (fPoints[pid*4+1]+GND_IMG_OFFY1)/(GND_IMG_DY1+0.000001);
        if(ix<0 || ix>=GND_IMG_NX1 || iy<0 || iy>=GND_IMG_NY1)
        {
            // 超出 pGndImg1 栅格范围的点在 tmpLabel1 中对应位置标记为-1
            tmpLabel1[pid]=-1;
            continue;
        }

        // 在范围内的点在 tmpLabel1 中记录其在 pGndImg1 中的索引
        int iid=ix+iy*GND_IMG_NX1;
        tmpLabel1[pid]=iid;

        // 找落在 pGndImg1 栅格中的z最小的点
        if(pGndImg1[iid]>fPoints[pid*4+2])
        {
            pGndImg1[iid]=fPoints[pid*4+2];
        }

    }

    // 记录地面点的数量
    int pnum=0;
    // 遍历输入进来的每一个点
    for(int pid=0;pid<pointNum;pid++)
    {
        if(tmpLabel1[pid]>=0) // >= 0 表明其在栅格范围内
        {
            // 在 pGndImg1 栅格中的，比最低点高0.4m范围内的点被标记为地面种子点
            if(pGndImg1[tmpLabel1[pid]]+0.4>fPoints[pid*4+2])
            {
                pLabel[pid]=1;
                pnum++;
            }
        }
    }
    free(pGndImg1);
    free(tmpLabel1);

    // 遍历输入进来的每一个点
    for(int pid=0;pid<pointNum;pid++)
    {
        // 在前一次比较中被标记为地面种子点
        if(pLabel[pid]==1) 
        {
            // 若该点的z值大于1，标记为非地面种子点     ****重要参数，z值高于1m的都不认为是地面点***
            if(fPoints[pid*4+2]>1)
            {
                pLabel[pid]=0;
            }
            // sqrt(x * x + y * y)小于100的，且z值大于0.5的也标记为非地面点
            else if(fPoints[pid*4]*fPoints[pid*4]+fPoints[pid*4+1]*fPoints[pid*4+1]<100)
            {
                if(fPoints[pid*4+2]>0.5)
                {
                    pLabel[pid]=0;
                }
            }

        }
        else
        {
            // sqrt(x * x + y * y)小于400的，且z值小于0.2的标记为地面点种子点
            if(fPoints[pid*4]*fPoints[pid*4]+fPoints[pid*4+1]*fPoints[pid*4+1]<400)
            {
                if(fPoints[pid*4+2]<0.2)
                {
                    pLabel[pid]=1;
                }
            }
        }

    }

    // 存20m范围内的地面点种子点，进行平面拟合
    pcl::PointCloud<pcl::PointXYZI>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZI>());
    int in_zone[pointNum] = {0};
    for(int pid=0;pid<pointNum;pid++)
    {
        // 只选20m范围内的点
        if(fPoints[pid*4]*fPoints[pid*4]+fPoints[pid*4+1]*fPoints[pid*4+1]<400)
        {
            in_zone[pid] = 1;
            if (pLabel[pid]==1)
            {
                pcl::PointXYZI p;
                p.x = fPoints[pid*4];
                p.y = fPoints[pid*4 + 1];
                p.z = fPoints[pid*4 + 2];
                p.intensity = fPoints[pid*4 + 3];
                cloud->points.push_back(p);
            }
        }
    }

    Eigen::Matrix3f cov;
	Eigen::Vector4f pc_mean;
    Eigen::MatrixXf normal_;

    // 通过计算地面点协方差矩阵的奇异值分解获取地面的法向量
	pcl::computeMeanAndCovarianceMatrix(*cloud, cov, pc_mean);
	Eigen::JacobiSVD<Eigen::MatrixXf> svd(cov,Eigen::DecompositionOptions::ComputeFullU);

    // 选择沿z轴的法向量
	normal_ = (svd.matrixU().col(2)); 
	Eigen::Vector3f seeds_mean = pc_mean.head<3>();

	//  normal.T * [x,y,z] = -d
    // 得到拟合平面的d
	float d_ = -(normal_.transpose()*seeds_mean)(0,0);

    // 获取拟合平面的额高度阈值
	float th_dist_d_ = 0.3 - d_;

    // 构造 n * 3 的点云矩阵，计算每一个点到平面的距离
    Eigen::MatrixXf points(pointNum, 3);
    for(int k = 0; k < pointNum; k++)
    {
        points.row(k) << fPoints[k*4], fPoints[k*4+ 1], fPoints[k*4+ 2];
    }

    // ground plane model   若到拟合平面的距离小于阈值，且范围在20m以内，则为最终的地面点
    Eigen::VectorXf result = points * normal_;

    for (int k = 0; k < pointNum; k++)
    {
        if (!in_zone[k])
            continue;
        if (result[k] < th_dist_d_)
        {
            pLabel[k] = 1;
        }
        else
        {
            pLabel[k] = 0;
        }
        
    }

    // 统计地面点数量
    gnum=0;
    for(int pid=0;pid<pointNum;pid++)
    {
        if(pLabel[pid]==1)
        {
            gnum++;
        }
    }

    return gnum;
}

