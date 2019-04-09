#define NUM 80
#define EXCHANGE_VALUE(a,b,temp) do {temp=a;a=b;b=temp;}while(0)
#define PRINT_ARRAY_TO_CONSOLE(array,st,ed,i,str) \
    do {\
        printf("%s\n", str);\
        for(i=st; i<=ed; i++){\
            printf("%.20f\n",array[i]);\
        }\
    } while(0)
#define FILTER struct Filter*
#define FILTER_DATA_FILENAME "filter_data.dat"  // フィルタ係数を出力するファイル
#define GAIN_DATA_FILENAME "gain_data.dat"      // フィルタの振幅特性を出力するファイル
#include <stdio.h>
#include <math.h>
struct Filter{
    int N;                  // 次数
    int k;                  // 実行回数
    double a[NUM];          // フィルタ係数
    double nfp;             // 通過域端周波数
    double nfs;             // 阻止域端周波数
    double nafp;            // 通過域端角周波数
    double nafs;            // 阻止域端角周波数
    double nafe[NUM];       // 極値周波数(k番目)
    double nafe_f[NUM];     // 生成された極値周波数(k+1番目)
    double max_error_gain;  // 最大誤差
    double diff_nafe_th;    // 極値周波数の変動の閾値
    double diff_seek_th;    // 極値探査における周波数変動の閾値
};
void print_result(FILTER);                      // 結果の出力
void initialize(char**, FILTER);                // フィルタの初期化
void refresh(FILTER);                           // フィルタのk+1回目への更新
void exec_remez(FILTER);                        // remezアルゴリズム実行
void seek_extremum_value(FILTER);               // 極値探査
int is_need_to_repeat(FILTER);                  // 繰り返すか否か判定する
double filter_gain(int, double, FILTER);        // フィルタのゲイン
double desired_filter_gain(double, FILTER);     // フィルタの所望ゲイン
double diff_filter_gain(int, double, FILTER);   // 微小時間経過によるゲインの差分

// コマンドライン引数の内容
// argv[1] : [int]    フィルタ係数の数(実際のフィルタ次数は2倍+1個になる)
// argv[2] : [double] 通過域端周波数
// argv[3] : [double] 阻止域端周波数

int main(int argc, char *argv[]){
    int i;
    struct Filter f;
    initialize(argv, &f);
    do{
        refresh(&f);
        exec_remez(&f);
        seek_extremum_value(&f);
    }while(is_need_to_repeat(&f));
    print_result(&f);
    return 0;
}

void print_result(FILTER f){
    int i,j, point_num;
    double g[3], l, h[NUM*2+3];
    FILE *fp;
    fp = fopen(FILTER_DATA_FILENAME,"wb");
    i = f->N*2+1;
    fwrite(&i, sizeof(unsigned int), 1, fp);
    for(j=0; j<f->N; j++){
        h[j] = f->a[f->N-j]/2.0;
    }
    h[f->N] = f->a[0];
    for(j=f->N+1; j<i; j++){
        h[j] = f->a[j-f->N]/2.0;
    }
    fwrite(h, sizeof(double), i, fp);
    fclose(fp);
    point_num = 100*f->N;
    l = M_PI/(double)point_num;
    i = point_num+1;
    fp = fopen(GAIN_DATA_FILENAME,"wb");
    fwrite(&i, sizeof(unsigned int), 1, fp);
    for(i = 0; i<=point_num; i++){
        g[0] = l*(double)i;
        g[1] = filter_gain(i, l, f);
        g[2] = 20*log10(fabs(g[1]));
        fwrite(g, sizeof(double), 3, fp);
    }
    fclose(fp);
    PRINT_ARRAY_TO_CONSOLE(h,0,f->N*2,i,"結果:フィルタ係数");
}

void initialize(char *argv[],FILTER f){
    int i, nafe_p_num, nafe_s_num;
    sscanf(argv[1],"%d",&f->N);
    sscanf(argv[2],"%lf",&f->nfp);
    sscanf(argv[3],"%lf",&f->nfs);
    f->k=-1;
    f->nafp = f->nfp*2*M_PI;
    f->nafs = f->nfs*2*M_PI;
    f->diff_nafe_th = pow(10.0,-8.0);
    f->diff_seek_th = pow(10.0,-8.0);
    // 極値周波数の初期値を設定
    nafe_p_num = (int)((double)(f->N-2)*f->nafp/(f->nafp+(M_PI-f->nafs))+0.5);
    nafe_s_num = (int)((double)(f->N-2)*(M_PI-f->nafs)/(f->nafp+(M_PI-f->nafs))+0.5);
    if(nafe_p_num+nafe_s_num > f->N-2){
        nafe_p_num--;
    }
    f->nafe_f[0] = 0.0;
    for(i=1; i<=nafe_p_num; i++){
        f->nafe_f[i] = (f->nafp/(double)(nafe_p_num+1.0)*(double)i);
    }
    f->nafe_f[nafe_p_num+1] = f->nafp;
    f->nafe_f[nafe_p_num+2] = f->nafs;
    for(i=1; i<=nafe_s_num; i++){
        f->nafe_f[nafe_p_num+2+i] = ((M_PI-f->nafs)/(double)(nafe_s_num+1.0)*(double)i+f->nafs);
    }
    f->nafe_f[f->N+1] = M_PI;
}

void refresh(FILTER f){
    int i;
    // 算出した極値周波数を別の配列に移す
    for(i=0; i<=f->N+1; i++){
        f->nafe[i] = f->nafe_f[i];
    }
    f->k++;
}

void exec_remez(FILTER f){
    int i,j,i2,max_num;
    double mat_A[NUM][NUM*2], vec_b[NUM], temp;
    // 行列Aの生成
    for(j=0; j<=f->N; j++){
        for(i=0; i<=f->N+1; i++){
            mat_A[i][j] = cos(j*f->nafe[i]);
        }
    }
    for(i=0; i<=f->N+1; i++){
        mat_A[i][f->N+1] = (double)(1-2*(i%2));
    }
    // ベクトルbの生成
    for(i=0; i<=f->N+1; i++){
        vec_b[i] = desired_filter_gain(f->nafe[i] ,f);
    }
    // 掃き出し法による解の算出
    for(i=0; i<=f->N; i++){
        //ピボット選択
        max_num = i;
        for(j=i; j<=f->N+1; j++){
            if(mat_A[j][i] > mat_A[max_num][i]){
                max_num = j;
            }
        }
        if(max_num != i){
            for(j=i; j<=f->N+1; j++){
                EXCHANGE_VALUE(mat_A[i][j], mat_A[max_num][j], temp);
            }
            EXCHANGE_VALUE(vec_b[i], vec_b[max_num], temp);
        }
        //掃き出し
        temp = mat_A[i][i];
        for(j=i; j<=f->N+1; j++){
            mat_A[i][j] /= temp;
        }
        vec_b[i] /= temp;
        for(i2=i+1; i2<=f->N+1; i2++){
            temp = mat_A[i2][i];
            for(j=i; j<=f->N+1; j++){
                mat_A[i2][j] -= mat_A[i][j]*temp;
            }
            vec_b[i2] -= vec_b[i]*temp;
        }
    }
    temp = mat_A[f->N+1][f->N+1];
    mat_A[f->N+1][f->N+1] /= temp;
    vec_b[f->N+1] /= temp;
    // 後退入力
    for(j=f->N+1; j>=1; j--){
        for(i=j-1; i>=0; i--){
            temp = mat_A[i][j]/mat_A[j][j];
            mat_A[i][j] -= mat_A[j][j]*temp;
            vec_b[i] -= vec_b[j]*temp;
        }
    }
    // 解の代入
    for(i=0; i<=f->N; i++){
        f->a[i] = vec_b[i];
    }
    f->max_error_gain = vec_b[f->N+1];
}

void seek_extremum_value(FILTER f){
    int i, j, n=0, point_num;
    double temp, nafm;
    double l, dfg, dfg_p;
    double naf_init[NUM]; 
    nafm = (f->nafp + f->nafs)/2.0;
    // 極値の探査
    point_num = 100*f->N;
    l = M_PI/(double)point_num;
    dfg_p = diff_filter_gain(1, l, f);
    for(i=2; i<point_num; i++){
        dfg = diff_filter_gain(i, l, f);
        if((dfg>0)-(dfg_p>0)){
            f->nafe_f[n] = l*(double)i;
            n++;
        }
        dfg_p = dfg;
    }
    f->nafe_f[n] = 0.0;
    f->nafe_f[n+1] = M_PI;
    f->nafe_f[n+2] = f->nafp;
    f->nafe_f[n+3] = f->nafs;
    n = n+3;
    // 昇順ソート
    for(i=0; i<=n; i++){
        for(j=i+1; j<=n; j++){
            if(f->nafe_f[j] < f->nafe_f[i]){
                EXCHANGE_VALUE(f->nafe_f[i], f->nafe_f[j], temp);
            }
        }
    }
}

int is_need_to_repeat(FILTER f){
    int n;
    double diff_nafe, diff_nafe_max = 0;
    // スレッショルド以下であれば終了(0を返す)
    for(n=0; n<=f->N+1; n++){
        diff_nafe = fabs(f->nafe_f[n] - f->nafe[n]);
        if(diff_nafe > diff_nafe_max){diff_nafe_max = diff_nafe;}
    }
    if(diff_nafe_max < f->diff_nafe_th){
        return 0;
    }else{
        return 1;
    }
}

double filter_gain(int point, double l, FILTER f){
    int n;
    double res = 0.0;
    for(n=f->N+1; n>=0; n--){
        res += f->a[n]*cos(fmod((double)n*((double)(point)*l),2.0*M_PI));
    }
    return res;
}

double diff_filter_gain(int point, double l, FILTER f){
    int n;
    double res = 0.0;
    for(n=f->N+1; n>=0; n--){
        res += f->a[n]*cos(fmod((double)n*((double)(point+1)*l),2.0*M_PI));
    }
    for(n=f->N+1; n>=0; n--){
        res -= f->a[n]*cos(fmod((double)n*((double)(point)*l),2.0*M_PI));
    }
    return res;
}

double desired_filter_gain(double naf, FILTER f){
    if (naf <= f->nafp){
        return 1.0;
    }else if (naf >= f->nafs){
        return 0.0;
    }else{
        return -1.0;
    }
}