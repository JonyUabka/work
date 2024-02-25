/**
*	@file	ISL_FileOperation.cpp
*	@brief	[Source file of File Operation Functions], 文件操作；
*	@see	ISeisLib Manual
*	@author [Liu Baihong, Yang Qiang, Song ZhiXiang], 刘百红、杨强、宋志翔；
*	@date	2014-07-02
*	@refer
*/

#include "ISL_FileOperation.h"


namespace ISLib {


/* 从二进制文件中读取浮点型数组 */
float * ISL_readFloatFromBinaryFile(char * fname, long num)
{
	if (fname == NULL)
		return NULL;

	float * val = new float[num];
	FILE * out;
	out = fopen(fname, "rb+");
	if (out == NULL) {
		printf("open file error !");
		return NULL;
	}

	for (long i = 0; i < num; i++) {
		fread(&val[i], sizeof(float), 1, out);
	}
	fclose(out);

	return val;
}

/* 将一个浮点型数组保存为二进制文件 */
int ISL_writeFloatToBinaryFile(char * fname, float * data, long num)
{
	if (fname == NULL || num <= 0 || data == NULL)
		return -1;

	ofstream out(fname, ios::binary);
	for (long i = 0; i < num; i++)
		out.write((char*) &data[i], sizeof(float));

	out.close();

	return 1;
}


/* 清空一个二进制文件 */
void ISL_clearBinaryFile(char * fname)
{
	FILE *fp = fopen(fname, "wb");
	if (fp == NULL)
		return;

	float * val = NULL;
	fwrite(val, 0, 0, fp);

	fclose(fp);
}



/* 从一个文本文件逐行读取数据并放入浮点数组中 */
float * ISL_readFloatFromTextFile(char * fname, long &num)
{
	if (fname == NULL) {
		num = 0;
		return NULL;
	}

	vector<float> data;
	data.clear();

	FILE * fp;
	fp = fopen(fname, "r");
	fseek(fp, 0, SEEK_SET);

	float v = 0;
	while (!feof(fp)) //打开源文件
	{
		v = 0;
		fscanf(fp, "%f", &v);
		data.push_back(v);
	}
	fclose(fp);

	num = data.size();
	float * val = new float[num];

	for (long i = 0; i < num; i++) {
		val[i] = data[i];
	}
	data.clear();

	return val;
}





/* 将一个浮点数组逐行写入一个文本文件中 */
void ISL_writeFloatToTextFile(char * fname, float * data, long num)
{
	if (fname == NULL || data == NULL || num <= 0)
		return;

	FILE * fp;
	fp = fopen(fname, "wt");

	for (long i = 0; i < num; i++) {
		fprintf(fp, "%.5f\n", data[i]);
	}
	fclose(fp);
}

/* 读取文件的大小 */
long ISL_getFileSize(FILE * f)
{
	long pos = ftell(f);
	fseek(f, 0, SEEK_END);
	long len = ftell(f);
	fseek(f, pos, SEEK_SET);
	return len;
}

}/* End Of namespace ISLIB */
