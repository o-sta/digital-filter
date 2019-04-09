#define PRINT_ARRAY_TO_CONSOLE(array,st,ed,i,str) \
    do {\
        printf("%s\n", str);\
        for(i=st; i<=ed; i++){\
            printf("%.20f\n",array[i]);\
        }\
    } while(0)
#define ORDER 500                   // フィルタ次数の最大数
#define DATA_SIZE 2500              // メモリに確保する音声信号の数
#define FILTER struct Filter*       // フィルタの構造体
#define SIGNAL struct Signal*       // 信号の構造体
#define WAVHDR struct WaveHeader    // wavファイルのヘッダを格納する構造体
#include <math.h>
#include <stdio.h>
struct Filter{
    unsigned int order; // フィルタ次数
    double a[ORDER];    // フィルタ係数
    double nfp1;        // 阻止域開始正規化周波数
    double nfp2;        // 阻止域終了正規化周波数
    double nafp1;       // 阻止域開始正規化角周波数
    double nafp2;       // 阻止域終了正規化角周波数
    FILE *fp;           // フィルタ係数を保存するファイルポインタ
};
struct WaveHeader{
    char riff[4];
    unsigned int chunk_size;
    char format[4];
    char fmt[4];
    unsigned int chunk_byte_size;
    unsigned short int audio_format;
    unsigned short int channel;
    unsigned int sampling;
    unsigned int bytes_per_second;
    unsigned short int block_size;
    unsigned short int bits_per_sample;
    char data[4];
    unsigned int bytes_of_data;
};
struct Signal{
    WAVHDR h;   // ヘッダファイル
    FILE *fp;   // ファイルポインタ
};
void exec_filter(SIGNAL, SIGNAL, FILTER);               // フィルタ実行
void initialize(int, char**, SIGNAL, SIGNAL, FILTER);   // 初期化
void dispose(SIGNAL, SIGNAL, FILTER);                   // ファイルクローズ
void setup_filter(FILTER);                              // フィルタの設定
void print_result(FILTER);                              // 結果の出力

// コマンドライン引数の内容
// argv[1] : [文字列]  入力音声ファイル(*.wav)
// argv[2] : [文字列]  出力(ノイズ除去後)音声ファイル(*.wav)
// argv[3] : [int]    フィルタ次数(必ず奇数を入力する)
// argv[4] : [double] 阻止域端周波数(最小値)
// argv[5] : [double] 阻止域端周波数(最大値)

int main(int argc, char *argv[]){
    struct Signal si, so;   // 入力・出力信号の構造体
    struct Filter f;        // フィルタの構造体
    initialize(argc, argv, &si, &so, &f);
    exec_filter(&si, &so, &f);
    dispose(&si, &so, &f);
    print_result(&f);
    return 0;
}

void print_result(FILTER f){
    int i;
    FILE *fp;
    fp = fopen("filter_BEF.dat","wb");
    i = f->order;
    fwrite(&i, sizeof(unsigned int), 1, fp);
    fwrite(f->a, sizeof(double), i, fp);
    fclose(fp);
    PRINT_ARRAY_TO_CONSOLE(f->a,0,f->order-1,i,"出力結果:フィルタ係数");
}

void initialize(int argc, char *argv[], SIGNAL si, SIGNAL so, FILTER f){
    si->fp = fopen(argv[1], "rb");
    so->fp = fopen(argv[2], "wb");
    fread(&si->h, sizeof(WAVHDR), 1, si->fp);
    fwrite(&si->h, sizeof(WAVHDR), 1, so->fp);
    sscanf(argv[3],"%d",&f->order);
    sscanf(argv[4],"%lf",&f->nfp1);
    sscanf(argv[5],"%lf",&f->nfp2);
    f->nafp1 = 2.0*M_PI*f->nfp1;
    f->nafp2 = 2.0*M_PI*f->nfp2;
    setup_filter(f);
}

void dispose(SIGNAL si, SIGNAL so, FILTER f){
    fclose(si->fp);
    fclose(so->fp);
}

void exec_filter(SIGNAL si, SIGNAL so, FILTER f){
    int i, j=0, k, data_num;
    short int data_in[DATA_SIZE];
    short int buffer[ORDER];
    short int data_out[DATA_SIZE]={0};
    double data;
    // 1. 入力ファイルから信号を <DATA_SIZE> 個取り出してメモリに保存する
    // 2. メモリ内の信号をフィルタに通す(フィルタの状態はbufferに格納される)
    // 3. 信号を全て出力ファイルに書き込む
    // 4. 1に戻り入力ファイルの信号を全て読み終えるまで繰り返す
    while (feof(si->fp)==0){
        data_num = (int)fread(data_in, sizeof(short int), DATA_SIZE, si->fp);
        for(i=0; i<data_num; i++){
            buffer[j] = data_in[i];
            data = 0.0;
            for(k=0; k<f->order; k++){
                data += f->a[k]*(double)buffer[(k+j)%f->order];
            }
            data_out[i] = (short int)data;
            j = (j+1)%f->order;
        }
        fwrite(data_out, sizeof(short int), (size_t)data_num, so->fp);
    }
}

void setup_filter(FILTER f){
    int i, transition, n;
    // 窓関数法によるBEFの係数設定
    transition = (f->order-1)/2;
    for(i=0; i<f->order; i++){
        n = i-transition;
        f->a[i] = 0.0;
        if(n == 0){
            f->a[i] = -2.0*f->nfp2 + 2.0*f->nfp1 + 1;
        }else{
            f->a[i] = -2.0*f->nfp2*sin(f->nafp2*(double)(n))/(f->nafp2*(double)(n))+2.0*f->nfp1*sin(f->nafp1*(double)(n))/(f->nafp1*(double)(n));
        }
    }
    // ハニング窓関数の重畳
    for(i=0; i<f->order; i++){
        n = i-transition;
        f->a[i] *= (0.5+0.5*cos((2.0*M_PI*n)/(f->order-1)));
    }
}