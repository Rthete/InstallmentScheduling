#ifndef _EXP_4_H
#define _EXP_4_H

#include "method.h"
#include <algorithm>
#include <fstream>
#include <sstream>

namespace exp_4 {
// 第一部分：对比三个模型MISRRLL, APMISRR, TolerMIS

// 不含容错_15台处理机_时间T-趟数m
void without_error_15_m();
// 不含容错_15台处理机_利用率ur-趟数m
void ur_without_error_15_m();
// 不含容错_15台处理机_时间T-任务量W
void without_error_15_W();
// 不含容错_15台处理机_利用率ur-任务量W
void ur_without_error_15_W();
// 不含容错_15台处理机_时间T-回传比例theta
void without_error_15_theta();
// 不含容错_15台处理机_利用率ur-回传比例theta
void ur_without_error_15_theta();

// 不含容错_30台处理机_时间T-趟数m
void without_error_30_m();
// 不含容错_30台处理机_利用率ur-趟数m
void ur_without_error_30_m();
// 不含容错_30台处理机_时间T-任务量W
void without_error_30_W();
// 不含容错_30台处理机_利用率ur-任务量W
void ur_without_error_30_W();
// 不含容错_30台处理机_时间T-回传比例theta
void without_error_30_theta();
// 不含容错_15台处理机_利用率ur-回传比例theta
void ur_without_error_30_theta();

// 第二部分：对比MISRRLL单模型

// 不含容错_15台处理机_图例任务量W_时间T-趟数m
void W_vs_m_15();
// 不含容错_15台处理机_图例回传比例theta_时间T-趟数m
void theta_vs_m_15();
// 不含容错_15台处理机_图例lambda1_时间T-趟数m
void lambda1_vs_m_15();
// 不含容错_15台处理机_图例lambda2_时间T-趟数m
void lambda2_vs_m_15();

// 不含容错_30台处理机_图例任务量W_时间T-趟数m
void W_vs_m_30();
// 不含容错_30台处理机_图例回传比例theta_时间T-趟数m
void theta_vs_m_30();
// 不含容错_30台处理机_图例lambda1_时间T-趟数m
void lambda1_vs_m_30();
// 不含容错_30台处理机_图例lambda2_时间T-趟数m
void lambda2_vs_m_30();

// 第三部分：对比MISRRLL单模型 lambda

// 不含容错_15台处理机_图例lambda1_时间T-趟数m
void labmda1_T_vs_m_15();
// 不含容错_15台处理机_图例lambda2_时间T-趟数m
void labmda2_T_vs_m_15();
// 不含容错_30台处理机_图例lambda1_时间T-趟数m
void labmda1_T_vs_m_30();
// 不含容错_30台处理机_图例lambda2_时间T-趟数m
void labmda2_T_vs_m_30();

// 第四部分：容错实验
void error_MISRRLL_30();
void ur_error_MISRRLL_30();
void error_cmp_3_models_ur();

// 第五部分：恢复实验
void recover_MISRRLL_15();
void recover_TolerMIS_15();
void recover_myAPMISRR_15();
void recover_cmp_3_models_ur();
void recover_cmp_3_models_ur_30();

} // namespace exp_4

#endif //_EXP_4_H