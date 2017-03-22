#include "stdafx.h"

#include "Parameterization.h"
#include <assert.h>
#include <math.h>
#include <mkl.h>
#include <algorithm>
#include <fstream>

CParameterization::CParameterization(CParameterization& other)
{

}

CParameterization::~CParameterization()
{
	/*
	for (int i = 0; i < m_innerVerticesResult.size(); i++)
	{
		delete m_innerVerticesResult[i];
	}
	for (int i = 0; i < m_borderVerticesResult.size(); i++)
	{
		delete m_borderVerticesResult[i];
	}
	for (int i = 0; i < m_nVertices; i++)
	{
		delete m_pVerticesTrans[i];
	}
	delete[] m_pVerticesTrans;
	*/
}

void CParameterization::BorderVerticesParameterize(int TYPE, float args)
{
	assert(TYPE == TYPE_SQUARE || TYPE == TYPE_CIRCLE);
	switch (TYPE)
	{
		case TYPE_SQUARE:
		{
			m_fBorderLength = 0.0f;
			std::vector<float> weight;
			weight.push_back(0);
			for (int i = 1; i < m_borderVertices.size(); i++)
			{
				float distance = 
					pow(m_pVertices[m_borderVertices[i]].m_vPosition.x - m_pVertices[m_borderVertices[i - 1]].m_vPosition.x, 2) +
					pow(m_pVertices[m_borderVertices[i]].m_vPosition.y - m_pVertices[m_borderVertices[i - 1]].m_vPosition.y, 2) +
					pow(m_pVertices[m_borderVertices[i]].m_vPosition.z - m_pVertices[m_borderVertices[i - 1]].m_vPosition.z, 2);
				distance = sqrtf(distance);
				m_fBorderLength += distance;
				weight.push_back(m_fBorderLength);
			}
			int end = m_borderVertices.size() - 1;
			float distance =
				pow(m_pVertices[m_borderVertices[end]].m_vPosition.x - m_pVertices[m_borderVertices[0]].m_vPosition.x, 2) +
				pow(m_pVertices[m_borderVertices[end]].m_vPosition.y - m_pVertices[m_borderVertices[0]].m_vPosition.y, 2) +
				pow(m_pVertices[m_borderVertices[end]].m_vPosition.z - m_pVertices[m_borderVertices[0]].m_vPosition.z, 2);
			distance = sqrtf(distance);
			m_fBorderLength += distance;
			for (int i = 0; i < m_borderVertices.size(); i++)
			{
				weight[i] = weight[i] / m_fBorderLength;
			}
			float myu = 0.25;
			float offset = args / 2.0f;
			for (int i = 0; i < weight.size(); i++)
			{
				if (weight[i] < myu)
				{
					CPoint3D* point = new CPoint3D(weight[i] * args * 4.0f - offset, 0.0f - offset, 0.0f);
					m_borderVerticesResult.push_back(point);
				}
				else if (weight[i] < myu * 2)
				{
					CPoint3D* point = new CPoint3D(args - offset, (weight[i] - myu) * args * 4.0f - offset, 0.0f);
					m_borderVerticesResult.push_back(point);
				}
				else if (weight[i] < myu * 3)
				{
					CPoint3D* point = new CPoint3D((3 * myu - weight[i]) * args * 4.0f - offset, args - offset, 0.0f);
					m_borderVerticesResult.push_back(point);
				}
				else if (weight[i] < myu * 4)
				{
					CPoint3D* point = new CPoint3D(0.0f - offset, (4 * myu - weight[i]) * args * 4.0f - offset, 0.0f);
					m_borderVerticesResult.push_back(point);
				}
			}
			break;
		}
		case TYPE_CIRCLE:
		{
			m_fBorderLength = 0.0f;
			std::vector<float> weight;
			weight.push_back(0);
			for (int i = 1; i < m_borderVertices.size(); i++)
			{
				float distance =
					pow(m_pVertices[m_borderVertices[i]].m_vPosition.x - m_pVertices[m_borderVertices[i - 1]].m_vPosition.x, 2) +
					pow(m_pVertices[m_borderVertices[i]].m_vPosition.y - m_pVertices[m_borderVertices[i - 1]].m_vPosition.y, 2) +
					pow(m_pVertices[m_borderVertices[i]].m_vPosition.z - m_pVertices[m_borderVertices[i - 1]].m_vPosition.z, 2);
				distance = sqrtf(distance);
				m_fBorderLength += distance;
				weight.push_back(m_fBorderLength);
			}
			int end = m_borderVertices.size() - 1;
			float distance =
				pow(m_pVertices[m_borderVertices[end]].m_vPosition.x - m_pVertices[m_borderVertices[0]].m_vPosition.x, 2) +
				pow(m_pVertices[m_borderVertices[end]].m_vPosition.y - m_pVertices[m_borderVertices[0]].m_vPosition.y, 2) +
				pow(m_pVertices[m_borderVertices[end]].m_vPosition.z - m_pVertices[m_borderVertices[0]].m_vPosition.z, 2);
			distance = sqrtf(distance);
			m_fBorderLength += distance;
			for (int i = 0; i < m_borderVertices.size(); i++)
			{
				weight[i] = weight[i] / m_fBorderLength;
			}
			float radius = args * 0.5f;
			float circle = 2.0f * PI;
			for (int i = 0; i < weight.size(); i++)
			{
				CPoint3D* point = new CPoint3D(radius * cosf(circle * weight[i]), radius * sinf(circle * weight[i]), 0.0f);
				m_borderVerticesResult.push_back(point);
			}
			break;
		}
	}
}


/*
void CParameterization::BorderVerticesParameterize(int TYPE, float args)
{
	assert(TYPE == TYPE_SQUARE || TYPE == TYPE_CIRCLE);
	switch (TYPE)
	{
	case TYPE_SQUARE:
	{
		int num = m_borderVertices.size() / 4;
		float gap = args / (float)num;
		float offset = args / 2.0f;
		for (int i = 0; i < num; i++)
		{
			CPoint3D* point = new CPoint3D(gap * i - offset, 0.0f - offset, 0.0f);
			m_borderVerticesResult.push_back(point);
		}
		for (int i = 0; i < num; i++)
		{
			CPoint3D* point = new CPoint3D(args - offset, gap * i - offset, 0.0f);
			m_borderVerticesResult.push_back(point);
		}
		for (int i = 0; i < num; i++)
		{
			CPoint3D* point = new CPoint3D(args - gap * i - offset, args - offset, 0.0f);
			m_borderVerticesResult.push_back(point);
		}
		for (int i = 0; i < m_borderVertices.size() - 3 * num; i++)
		{
			CPoint3D* point = new CPoint3D(0.0f - offset, args - gap * i - offset, 0.0f);
			m_borderVerticesResult.push_back(point);
		}
		break;
	}
	case TYPE_CIRCLE:
		assert(false);
		break;
	}
}
*/

void CParameterization::InnerVerticesParameterize(int TYPE)
{
	assert(TYPE == TYPE_AVERAGE_PARAMETERIZATION
		|| TYPE == TYPE_SHAPE_PRESERVING_PARAMETERIZATION);

	int nInenerNumber = m_innerVertices.size();
	int nBorderNumber = m_borderVertices.size();
	int nTotalNumber = nInenerNumber + nBorderNumber;

	//au = b	a(N*N)  u(N*2) b(N*2) --> u = Inverse(a) * b
	
	float* fpMatrixA = new float[nTotalNumber * nTotalNumber];
	float* fpMatrixB = new float[nTotalNumber * 2];

	memset(fpMatrixA, 0, sizeof(float) * nTotalNumber * nTotalNumber);
	memset(fpMatrixB, 0, sizeof(float) * nTotalNumber * 2);

	//初始化Ai,i = 1
	for (int i = 0; i < nTotalNumber; i++)
	{ 
		setMatrix(fpMatrixA, nTotalNumber, i, i, 1.0f); 
	}

	//对于每一个内部的点，检索它周围的点
	for (int i = 0; i < nInenerNumber; i++)
	{
		int nNeighborsNumber = m_pVertices[m_innerVertices[i]].m_nNeighbor;
		//获取这个点的保形参数（只含有非零元素）
		std::vector<float>* myus;
		switch (TYPE)
		{
			case TYPE_AVERAGE_PARAMETERIZATION:
			{
				myus = new std::vector<float>;
				for (int w = 0; w < nNeighborsNumber; w++)
				{
					myus->push_back(1.0f / nNeighborsNumber);
				}
				break;
			}

			case TYPE_SHAPE_PRESERVING_PARAMETERIZATION:
			{
				myus = ShapePreservingParameterization(m_innerVertices[i]);
				break;
			}

			default:
				assert(false);
				break;
		}

		//ADD
		//计算（9）式等号右边矩阵的每一行（0-i）
		CPoint3D b = CPoint3D(0.0f, 0.0f, 0.0f);
		
		for (int j = 0; j < nNeighborsNumber; j++)
		{
			//找到邻点对应的源数据中的索引号
			int nNeighberIndex = m_pVertices[m_innerVertices[i]].m_npNeighborVertexIndices[j];
			//计算在经过分组（内部，外部点）后，邻点对应的序号1-N
			int nInsertedColumn;
			std::vector<int>::iterator iter;
			iter = std::find(m_innerVertices.begin(), m_innerVertices.end(), nNeighberIndex);
			//内部点在前
			if (iter != m_innerVertices.end())
			{
				nInsertedColumn = std::distance(m_innerVertices.begin(), iter);
			}
			else
			{
				//内部点没有，外部点中查找
				iter = std::find(m_borderVertices.begin(), m_borderVertices.end(), nNeighberIndex);
				if (iter != m_borderVertices.end())
				{
					nInsertedColumn = std::distance(m_borderVertices.begin(), iter)
						+ m_innerVertices.size();
				}
				else
				{
					printf("------------------------------------------------------------------");
					printf("error occur during the process of inner vertices parameterization:");
					printf("could not find the adjacent vertex!");
					printf("------------------------------------------------------------------");
					delete[] fpMatrixA;
					delete[] fpMatrixB;
					assert(0);
				}
			}

			//ADD
			//如果邻点为外部点，累加（9）右侧
			int index = nInsertedColumn - m_innerVertices.size();
			if (index >= 0){
				b.x += (*myus)[j] * m_borderVerticesResult[index]->x;
				b.y += (*myus)[j] * m_borderVerticesResult[index]->y;
			}

			setMatrix(fpMatrixA, nTotalNumber, i, nInsertedColumn, -(*myus)[j]);
		}
		delete myus;
		//ADD
		//在这个点处理完成后（9）等号右侧的累加也完成，将它添加到矩阵B中
		setMatrix(fpMatrixB, 2, i, 0, b.x);
		setMatrix(fpMatrixB, 2, i, 1, b.y);
	}

	//解方程
	int* ipiv = new int[nTotalNumber * nTotalNumber];
	int result = LAPACKE_sgesv(LAPACK_ROW_MAJOR, nTotalNumber, 2, fpMatrixA, nTotalNumber, ipiv, fpMatrixB, 2);
	
	if (result == 0)
	{
		for (int i = 0; i < nInenerNumber; i++)
		{
			CPoint3D* solution = new CPoint3D(getMatrix(fpMatrixB, 2, i, 0), getMatrix(fpMatrixB,2 ,i, 1), 0.0f);
			m_innerVerticesResult.push_back(solution);
		}
	}
	else
	{
		printf("------------------------------------------------------------------");
		printf("error occur during the process of solving the linear equation:");
		printf("error code:%d", result);
		printf("If info = -i, parameter i had an illegal value.");
		printf("If info = i,U(i,i)(computed in double precision for mixed precision");
		printf("subroutines) is exactly zero. The factorization has been completed,");
		printf("but the factor U is exactly singular, so the solution could not be ");
		printf("computed.");
		printf("------------------------------------------------------------------");
		delete[] ipiv;
		delete[] fpMatrixA;
		delete[] fpMatrixB;
		assert(0);
	}

	delete[] ipiv;
	delete[] fpMatrixA;
	delete[] fpMatrixB;
}

/*
void CParameterization::DataInitialize()
{
	bool* bpEdgeTag = new bool[m_nEdges];
	for (int i = 0; i < m_nEdges; i++)
	{
		bpEdgeTag[i] = false;
	}

	//点
	for (int i = 0; i < m_nVertices; i++)
	{
		HE_vert* tmp_vert = new HE_vert();
		tmp_vert->x = m_pVertices[i].m_vPosition.x;
		tmp_vert->y = m_pVertices[i].m_vPosition.y;
		tmp_vert->z = m_pVertices[i].m_vPosition.z;
		vertices.push_back(tmp_vert);
	}

	//面和边
	for (int i = 0; i < m_nTriangles; i++)
	{
		HE_face* tmp_face = new HE_face();
		faces.push_back(tmp_face);
		HE_edge* last = NULL;
		HE_edge* first = NULL;
		for (int j = 0; j < 3; j++)
		{
			HE_vert* tmp_vert = NULL;
			if (!bpEdgeTag[i * 3 + j])
			{
				tmp_vert = vertices[m_pEdges[m_pTriangles[i].m_nEdgeIndices[j]].m_nVertexIndices[0]];
				bpEdgeTag[i * 3 + j] = true;
			}
			else
			{
				tmp_vert = vertices[m_pEdges[m_pTriangles[i].m_nEdgeIndices[j]].m_nVertexIndices[1]];
			}
			HE_edge* tmp_edge = new HE_edge();
			tmp_edge->face = tmp_face;
			tmp_edge->vert = tmp_vert;
			tmp_vert->edge = tmp_edge;
			if (first == NULL)
			{
				first = tmp_edge;
				tmp_face->edge = tmp_edge;
			}
			if (last != NULL)
			{
				last->next = tmp_edge;
			}
			last = tmp_edge;
			edges.push_back(tmp_edge);
		}
		last->next = first;
	}

	for (int i = 0; i < m_nTriangles; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			HE_face* other_face = faces[m_pEdges[m_pTriangles[i].m_nEdgeIndices[j]].m_nTriangleIndices[1]];
			HE_face* this_face = faces[i];
			
		}
	}

	delete[] bpEdgeTag;
}
*/

void CParameterization::FindBorderInnerVertices()
{
	bool* tag = new bool[m_nVertices];
	for (int i = 0; i < m_nVertices; i++)
	{
		tag[i] = false;
	}
	//先找到边界点
	std::vector<int> borderVertices;
	for (int i = 0; i < m_nEdges; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			if (m_pEdges[i].m_nTriangleIndices[0] == -1 || m_pEdges[i].m_nTriangleIndices[1] == -1)
			{
				if (!tag[m_pEdges[i].m_nVertexIndices[0]])
				{
					borderVertices.push_back(m_pEdges[i].m_nVertexIndices[0]);
					tag[m_pEdges[i].m_nVertexIndices[0]] = true;
				}
				if (!tag[m_pEdges[i].m_nVertexIndices[1]])
				{
					borderVertices.push_back(m_pEdges[i].m_nVertexIndices[1]);
					tag[m_pEdges[i].m_nVertexIndices[1]] = true;
				}
			}
		}
	}
	delete[] tag;

	//边界点按顺序排好
	tag = new bool[m_nEdges];
	for (int i = 0; i < m_nEdges; i++)
	{
		tag[i] = false;
	}
	int nFind = 1, cur = borderVertices[0];
	while (nFind <= borderVertices.size())
	{
		int index = cur;
		m_borderVertices.push_back(index);
		for (int i = 0; i < m_nEdges; i++)
		{
			if (!tag[i])
			{
				if (m_pEdges[i].m_nVertexIndices[0] == index || m_pEdges[i].m_nVertexIndices[1] == index)
				{
					if (m_pEdges[i].m_nTriangleIndices[0] == -1 || m_pEdges[i].m_nTriangleIndices[1] == -1)
					{
						if (m_pEdges[i].m_nVertexIndices[0] == index)
						{
							cur = m_pEdges[i].m_nVertexIndices[1];
						}
						else
						{
							cur = m_pEdges[i].m_nVertexIndices[0];
						}
						tag[i] = true;
						nFind++;
						break;
					}
				}
			}
		}
	}
	delete[] tag;

	for (int i = 0; i < m_nVertices; i++)
	{
		if (std::find(m_borderVertices.begin(), m_borderVertices.end(), i) == m_borderVertices.end())
		{
			m_innerVertices.push_back(i);
		}
	}
}

//std::vector<float>* CParameterization::ShapePreservingParameterization(int i)
//{
//	int nNeighborNumber = m_pVertices[i].m_nNeighbor;
//	float* fpAngles = new float[nNeighborNumber];
//	float fSumAngle = 0.0f;
//
//	//计算i点周围邻边夹角的权重（0-1）
//	fpAngles[0] = 0.0f;
//	for (int j = 1; j < nNeighborNumber; j++)
//	{
//		CPoint3D pa = CPoint3D(
//			m_pVertices[m_pVertices[i].m_npNeighborVertexIndices[j - 1]].m_vPosition.x - m_pVertices[i].m_vPosition.x,
//			m_pVertices[m_pVertices[i].m_npNeighborVertexIndices[j - 1]].m_vPosition.y - m_pVertices[i].m_vPosition.y,
//			m_pVertices[m_pVertices[i].m_npNeighborVertexIndices[j - 1]].m_vPosition.z - m_pVertices[i].m_vPosition.z
//			);
//		CPoint3D pb = CPoint3D(
//			m_pVertices[m_pVertices[i].m_npNeighborVertexIndices[j]].m_vPosition.x - m_pVertices[i].m_vPosition.x,
//			m_pVertices[m_pVertices[i].m_npNeighborVertexIndices[j]].m_vPosition.y - m_pVertices[i].m_vPosition.y,
//			m_pVertices[m_pVertices[i].m_npNeighborVertexIndices[j]].m_vPosition.z - m_pVertices[i].m_vPosition.z
//			);
//		float fAngleCos = (pa & pb) / (pa.length() * pb.length());
//		float fAngle = acosf(fAngleCos);	
//		fpAngles[j] = fAngle;
//		fSumAngle += fAngle;
//	}
//	CPoint3D pa = CPoint3D(
//		m_pVertices[m_pVertices[i].m_npNeighborVertexIndices[nNeighborNumber - 1]].m_vPosition.x - m_pVertices[i].m_vPosition.x,
//		m_pVertices[m_pVertices[i].m_npNeighborVertexIndices[nNeighborNumber - 1]].m_vPosition.y - m_pVertices[i].m_vPosition.y,
//		m_pVertices[m_pVertices[i].m_npNeighborVertexIndices[nNeighborNumber - 1]].m_vPosition.z - m_pVertices[i].m_vPosition.z
//		);
//	CPoint3D pb = CPoint3D(
//		m_pVertices[m_pVertices[i].m_npNeighborVertexIndices[0]].m_vPosition.x - m_pVertices[i].m_vPosition.x,
//		m_pVertices[m_pVertices[i].m_npNeighborVertexIndices[0]].m_vPosition.y - m_pVertices[i].m_vPosition.y,
//		m_pVertices[m_pVertices[i].m_npNeighborVertexIndices[0]].m_vPosition.z - m_pVertices[i].m_vPosition.z
//		);
//	float fAngleCos = (pa & pb) / (pa.length() * pb.length());
//	fSumAngle += acosf(fAngleCos);
//	
//	//计算转换后每一条边与x轴的夹角
//	for (int j = 0; j < nNeighborNumber; j++)
//	{
//		fpAngles[j] = 2 * PI * fpAngles[j] / fSumAngle;
//	}
//
//	//实现的时候通常假定 p = 0 and p1 = (||xjk − xi||, 0)
//	//先计算邻边的长度
//	float* fpLengths = new float[nNeighborNumber];
//	for (int j = 0; j < nNeighborNumber; j++)
//	{
//		CPoint3D v = CPoint3D(
//			m_pVertices[m_pVertices[i].m_npNeighborVertexIndices[j]].m_vPosition.x - m_pVertices[i].m_vPosition.x,
//			m_pVertices[m_pVertices[i].m_npNeighborVertexIndices[j]].m_vPosition.y - m_pVertices[i].m_vPosition.y,
//			m_pVertices[m_pVertices[i].m_npNeighborVertexIndices[j]].m_vPosition.z - m_pVertices[i].m_vPosition.z
//			);
//		fpLengths[j] = v.length();
//	}
//
//	//根据极坐标计算变换后的坐标
//	std::vector<CPoint3D*> transformedPoint;
//	for (int j = 0; j < nNeighborNumber; j++)
//	{
//		CPoint3D* point = new CPoint3D(fpLengths[j] * cosf(fpAngles[j]), fpLengths[j] * sinf(fpAngles[j]), 0.0f);
//		transformedPoint.push_back(point);
//	}
//
//	std::vector<float>* myus = new std::vector<float>[nNeighborNumber];
//
//	//实现的时候通常假定 p = 0
//	CPoint3D v1, v2, v3;
//	CPoint3D center = CPoint3D(0.0f, 0.0f, 0.0f);
//	int nIndexP1, nIndexP2, nIndexP3;
//	for (int j = 0; j < nNeighborNumber; j++)
//	{
//		//先算出中心点到邻边j的向量v3
//		int nCenterIndex = i;
//		v3 = (*transformedPoint[j]) - center;
//		nIndexP3 = j;
//		int nTagP1, nTagP2;
//		float fMinDelta = 9999.9f;
//		//重心坐标参数
//		float fU = 0.0f, fV = 0.0f, fW = 0.0f;
//		//搜索三角形的另外两个点
//		for (int k = 0; k < nNeighborNumber; k++)
//		{
//			//先计算中心点到俩点的向量v1v2
//			if (k == 0)
//			{
//				v1 = (*transformedPoint[nNeighborNumber - 1]) - center;
//				v2 = (*transformedPoint[0]) - center;
//				nIndexP1 = nNeighborNumber - 1;
//				nIndexP2 = 0;
//			}
//			else
//			{
//				v1 = (*transformedPoint[k - 1]) - center;
//				v2 = (*transformedPoint[k]) - center;
//				nIndexP1 = k - 1;
//				nIndexP2 = k;
//			}
//			if (nIndexP1 != nIndexP3 && nIndexP2 != nIndexP3)
//			{
//				//计算j点和找到的两个点组成三角形的面积
//				CPoint3D edge31 = (*transformedPoint[nIndexP1]) - (*transformedPoint[nIndexP3]);
//				CPoint3D edge21 = (*transformedPoint[nIndexP2]) - (*transformedPoint[nIndexP3]);
//				float fTriangleArea = 0.5 * (edge21 * edge31).length();		//S = 0.5*sin(aob) * |oa| * |ob|
//				//计算中心和三个点任意两点组成三角形的面积
//				float fAreaV1V2 = 0.5 * (v1 * v2).length(); 
//				float fAreaV1V3 = 0.5 * (v1 * v3).length();
//				float fAreaV2V3 = 0.5 * (v2 * v3).length();
//				float fArea = fAreaV1V2 + fAreaV1V3 + fAreaV2V3;
//				//如果三个三角形的面积更大则说明点在外侧
//				//尽量找在i点在三角内的
//				if (fArea - fTriangleArea < fMinDelta)
//				{
//					if (fArea != 0){
//						fMinDelta = fArea - fTriangleArea;
//						nTagP1 = nIndexP1;
//						nTagP2 = nIndexP2;
//						fU = fAreaV1V2 / fArea;
//						fV = fAreaV1V3 / fArea;
//						fW = fAreaV2V3 / fArea;
//					}
//				}
//			}
//		}
//		assert(fW >= 0 && fV >= 0 && fU >= 0);
//		assert(1.0f - fW - fV - fU < 0.001f && 1.0f - fW - fV - fU > -0.001f);
//		CPoint3D offset = (*transformedPoint[nIndexP3]) * fU + (*transformedPoint[nIndexP2]) * fV + (*transformedPoint[nIndexP1]) * fW;
//		myus[nTagP1].push_back(fW);
//		myus[nTagP2].push_back(fV);
//		myus[nIndexP3].push_back(fU);
//	}
//	//（15）式子
//	std::vector<float>* result = new std::vector<float>(nNeighborNumber, 0.0f);
//	for (int j = 0; j < nNeighborNumber; j++)
//	{
//		for (int k = 0; k < myus[j].size(); k++)
//		{
//			(*result)[j] += myus[j][k];
//		}
//		(*result)[j] = (*result)[j] / nNeighborNumber;
//	}
//
//	for (int i = 0; i < transformedPoint.size(); i++)
//	{
//		delete transformedPoint[i];
//	}
//	delete[] myus;
//	delete[] fpLengths;
//	delete[] fpAngles;
//	return result;
//}

std::vector<float>* CParameterization::ShapePreservingParameterization(int i)
{
	int nNeighborNumber = m_pVertices[i].m_nNeighbor;
	float* fpWeight = new float[nNeighborNumber];
	float fTotal = 0.0f;

	std::vector<float>* result = new std::vector<float>();
	for (int j = 0; j < nNeighborNumber; j++)
	{
		CPoint3D pji = CPoint3D(
			m_pVertices[m_pVertices[i].m_npNeighborVertexIndices[j]].m_vPosition.x - m_pVertices[i].m_vPosition.x,
			m_pVertices[m_pVertices[i].m_npNeighborVertexIndices[j]].m_vPosition.y - m_pVertices[i].m_vPosition.y,
			m_pVertices[m_pVertices[i].m_npNeighborVertexIndices[j]].m_vPosition.z - m_pVertices[i].m_vPosition.z
			);
		float fl = 1.0f / pji.length();
		
		fpWeight[j] = fl;
		fTotal += fl;
	}
	
	for (int j = 0; j < nNeighborNumber; j++)
	{
		result->push_back(fpWeight[j] / fTotal);
	}
	
	delete[] fpWeight;
	return result;
}

CPoint3D** CParameterization::Calculate(float args)
{
	FindBorderInnerVertices();
	BorderVerticesParameterize(TYPE_SQUARE, args);
	//BorderVerticesParameterize(TYPE_CIRCLE, args);
	InnerVerticesParameterize(TYPE_SHAPE_PRESERVING_PARAMETERIZATION);
	//InnerVerticesParameterize(TYPE_AVERAGE_PARAMETERIZATION);
	m_pVerticesTrans = new CPoint3D*[m_nVertices];
	for (int i = 0; i < m_innerVertices.size(); i++)
	{
		m_pVerticesTrans[m_innerVertices[i]] = new CPoint3D(
			m_innerVerticesResult[i]->x,
			m_innerVerticesResult[i]->y,
			m_innerVerticesResult[i]->z
			);
	}
	for (int i = 0; i < m_borderVertices.size(); i++)
	{
		m_pVerticesTrans[m_borderVertices[i]] = new CPoint3D(
			m_borderVerticesResult[i]->x,
			m_borderVerticesResult[i]->y,
			m_borderVerticesResult[i]->z
			);
	}
	return m_pVerticesTrans;
}

void CParameterization::Draw(float args)
{
	if (m_pVerticesTrans == nullptr)
	{
		Calculate(args);
		printf("calculated....\n");
	}
	if (m_pVerticesTrans != nullptr)
	{
		glDisable(GL_TEXTURE_2D);
		glDisable(GL_LIGHTING);
		glColor3f(0.5f, 0.5f, 0.5f);
		glBegin(GL_LINES);
		CPoint3D* v[2];
		for (int i = 0; i < m_nEdges; i++)
		{

			v[0] = m_pVerticesTrans[m_pEdges[i].m_nVertexIndices[0]];
			v[1] = m_pVerticesTrans[m_pEdges[i].m_nVertexIndices[1]];
			glVertex3f(v[0]->x, v[0]->y, v[0]->z);
			glVertex3f(v[1]->x, v[1]->y, v[1]->z);

		}
		glEnd();
	}
}

void* CParameterization::__pImage = NULL;
GLint CParameterization::__width = -1;
GLint CParameterization::__height = -1;

void CParameterization::Draw(float args, const char* filename)
{
	if (m_pVerticesTrans == nullptr)
	{
		Calculate(args);
		printf("calculated....\n");
	}

	loadImage(filename);
	
	//设置纹理过滤
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);


	//设置纹理环境
	//glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
	glEnable(GL_TEXTURE_2D);

	float offset = args / 2.0f;
	glBegin(GL_TRIANGLES);
	for (int i = 0; i < m_nTriangles; i++)
	{
		int index0 = m_pTriangles[i].m_nVertexIndices[0];
		int index1 = m_pTriangles[i].m_nVertexIndices[1];
		int index2 = m_pTriangles[i].m_nVertexIndices[2];

		CPoint3D v0 = CPoint3D(m_pVertices[index0].m_vPosition.x, m_pVertices[index0].m_vPosition.y, m_pVertices[index0].m_vPosition.z);
		CPoint3D v1 = CPoint3D(m_pVertices[index1].m_vPosition.x, m_pVertices[index1].m_vPosition.y, m_pVertices[index1].m_vPosition.z);
		CPoint3D v2 = CPoint3D(m_pVertices[index2].m_vPosition.x, m_pVertices[index2].m_vPosition.y, m_pVertices[index2].m_vPosition.z);
		
		CPoint3D _v0 = CPoint3D((m_pVerticesTrans[index0]->x + offset) / args, (m_pVerticesTrans[index0]->y + offset) / args, (m_pVerticesTrans[index0]->z + offset) / args);
		CPoint3D _v1 = CPoint3D((m_pVerticesTrans[index1]->x + offset) / args, (m_pVerticesTrans[index1]->y + offset) / args, (m_pVerticesTrans[index1]->z + offset) / args);
		CPoint3D _v2 = CPoint3D((m_pVerticesTrans[index2]->x + offset) / args, (m_pVerticesTrans[index2]->y + offset) / args, (m_pVerticesTrans[index2]->z + offset) / args);

		glTexCoord2f(_v0.x, _v0.y);  glVertex3f(v0.x, v0.y, v0.z);
		glTexCoord2f(_v1.x, _v1.y);  glVertex3f(v1.x, v1.y, v1.z);
		glTexCoord2f(_v2.x, _v2.y);  glVertex3f(v2.x, v2.y, v2.z);
	}
	glEnd();
}

void CParameterization::loadImage(const char *filename)
{
	//读取图像文件
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	if (__pImage == nullptr)
	{
		__pImage = gltloadbmp(filename, &__width, &__height);
	}
	if (__pImage)
	{
		//与glTexImage2D功能相似，都是生成纹理不同的是glTexImage2D只支持64*64、128*128、256*256的位图，而gluBuild2DMipmaps支持任意分辨率的位图   
		gluBuild2DMipmaps(GL_TEXTURE_2D, 3, __width, __height, GL_BGR_EXT, GL_UNSIGNED_BYTE, __pImage);
	}
}

GLbyte* CParameterization::gltloadbmp(const char *filename, GLint *width, GLint *height)
{
	//打开文件   
	FILE *pfile = fopen(filename, "rb");
	if (pfile == 0)
	{
		exit(0);
	}
	//读取图像大小 bmp图像宽高放在前两个字节？   
	fseek(pfile, 0x0012, SEEK_SET);//pfile指向偏移0x0012位，不成功则指向SEEK_SET即文件开头   
	fread(width, sizeof(*width), 1, pfile);//最多读1个元素，每个元素sizeof(*iwidth)个字节，iwidth用于接收数据的内存地址   
	fread(height, sizeof(*height), 1, pfile);
	//计算像素长度  
	GLint pixellength = (*width) * 3;//3:BGR（blue、green、red）   
	while (pixellength % 4 != 0)//读取的时候是按int读取的，也就是4字节   
	{
		pixellength++;
	}
	pixellength = pixellength*(*height);
	//读取像素数据  
	GLbyte *pixeldata = (GLbyte*)malloc(pixellength);
	if (pixeldata == 0)
	{
		exit(0);
	}
	fseek(pfile, 54, SEEK_SET);
	fread(pixeldata, pixellength, 1, pfile);
	//关闭文件  
	fclose(pfile);
	return pixeldata;
}

float CParameterization::getMatrix(float* matrix, int ld, int i, int j)
{
	assert(ld >= 1);
	assert(i >= 0 && j >= 0);
	return matrix[i * ld + j];
}

void CParameterization::setMatrix(float* matrix, int ld, int i, int j, float data)
{
	assert(ld >= 1);
	assert(i >= 0 && j >= 0);
	matrix[i * ld + j] = data;
}