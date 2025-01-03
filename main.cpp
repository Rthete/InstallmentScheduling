/*
 * @FilePath: \InstallmentScheduling\main.cpp
 * @Description:
 * @Author: rthete
 * @Date: 2023-03-14 16:19:26
 * @LastEditTime: 2024-01-25 14:39:13
 */

#include "MISRR.h"
#include "include/exp_1.h"
#include "include/exp_2.h"
#include "include/exp_3.h"
#include "include/exp_4.h"
#include "include/method.h"

int main() {
  /**
   * @brief 调度示例测试
   *
   */
  // ModelRunner::run_MISRR(4, 4, 400, 0.2, "../data/4-servers-example/");
  // ModelRunner::run_MISRR(4, 4, 400, 0.2, "../data/4-servers-example/", {3},
  // 2);

  /**
   * @brief 容错测试
   *
   */
  // ModelRunner::run_MISRR(30, 24, 8000, 0.3, "../data/exp1-30-servers/", {6,
  // 12}, 7); ModelRunner::run_MISRR(30, 24, 8000, 0.3,
  // "../data/exp1-30-servers/", {6, 12}, 23); ModelRunner::run_MISRR(30, 24,
  // 8000, 0.3, "../data/exp1-30-servers/"); ModelRunner::run_myAPMISRR(30, 0.5,
  // 24, 8000, 0.3, "../data/exp1-30-servers/", {6, 12}, 7);

  /**
   * @brief 带容错的三个模型对比实验
   *
   */
  // exp_3::error_SIS_15();
  // exp_3::error_TolerMIS_15();
  // exp_3::error_APMISRR_15();
  // exp_3::error_SIS_30();
  // exp_3::error_TolerMIS_30();
  // exp_3::error_APMISRR_30();
  // exp_3::error_TolerMIS_30_conflict();
  // exp_3::error_cmp_3_models_ur();
  // exp_3::error_save_each_time();
  // exp_3::error_TolerMIS_30_conflict_inst();
  // exp_3::error_TolerMIS_30_conflict_workload();

  /**
   * @brief 不带容错的单模型实验
   */

  // exp_2::W_vs_m_15();
  // exp_2::theta_vs_m_15();
  // exp_2::W_vs_m_30();
  // exp_2::theta_vs_m_30();

  /**
   * @brief 不带容错的三个模型实验
   */
  // without_error_15_m();
  // ur_without_error_15_m();
  // without_error_15_W();
  // ur_without_error_15_W();
  // without_error_15_theta();
  // ur_without_error_15_theta();

  // without_error_30_m();
  // ur_without_error_30_m();
  // without_error_30_W();
  // ur_without_error_30_W();
  // without_error_30_theta();
  // ur_without_error_30_theta();

  /**
   * @brief 不带容错的三个模型实验MISRRLL
   */
  // exp_4::without_error_15_m();
  // exp_4::ur_without_error_15_m();
  // exp_4::without_error_15_W();
  // exp_4::ur_without_error_15_W();
  // exp_4::without_error_15_theta();
  // exp_4::ur_without_error_15_theta();

  // exp_4::without_error_30_m();
  // exp_4::ur_without_error_30_m();
  // exp_4::without_error_30_W();
  // exp_4::ur_without_error_30_W();
  // exp_4::without_error_30_theta();
  // exp_4::ur_without_error_30_theta();

  /**
   * @brief 不带容错单模型实验MISRRLL
   */
  // exp_4::W_vs_m_15();
  // exp_4::W_vs_m_30();
  // exp_4::theta_vs_m_15();
  // exp_4::theta_vs_m_30();
  // exp_4::lambda1_vs_m_15();
  // exp_4::lambda1_vs_m_30();
  // exp_4::lambda2_vs_m_15();
  // exp_4::lambda2_vs_m_30();

  /**
   * @brief 不带容错单模型实验lambda MISRRLL
   */
  // exp_4::labmda1_T_vs_m_15();
  // exp_4::labmda1_T_vs_m_30();
  // exp_4::labmda2_T_vs_m_15();
  // exp_4::labmda2_T_vs_m_30();

  /**
   * @brief 带容错单模型实验MISRRLL
   */
  // exp_4::error_MISRRLL_30();
  // exp_4::ur_error_MISRRLL_30();
  // exp_4::error_cmp_3_models_ur();

  /**
   * @brief 带恢复单模型实验MISRRLL
   */
  // exp_4::recover_MISRRLL_15();
  // exp_4::recover_TolerMIS_15();
  exp_4::recover_myAPMISRR_15();

  // MISRRLL misrrll(15, 0.3, 24);
  // cout << misrrll.findFirstPositiveRealRoot(3, -6, -1200, 3000, -2000) <<
  // endl;

  // with_error_15();

  // ModelRunner::run_SIS();
  // ModelRunner::run_myAPMISRR(0.6, 8);
  // ModelRunner::run_MISRR(0.3, 8);
  // ModelRunner::run_PMIS();
  // ModelRunner::run_APMISRR(0.6);
  // ModelRunner::run_APMISRR_cost(0.6, 8);
  // ModelRunner::run_MISRRL(1, 8);
  // ModelRunner::run_MISRRLL(1, 0.7, 8);

  // ModelRunner::run_MISRRLL(1, 1, 34, 8000);
  // ModelRunner::run_MISRR(15, 34, 8000, 0.3, "../data/exp1-15-servers/");
  // ModelRunner::run_MISRR(15, 10, 5000, 0.3, "../data/15-servers-w-20/", {6,
  // 12}, 5); ModelRunner::run_MISRR(3, 3, 300, 0.3,
  // "../data/3-servers-example/");

  // test_MISRR_theta();
  // test_MISRR_all();
  // test_APMISRR_totalTime();
  // test_APMISRR_alpha();
  // test_APMISRR_beta();
  // test_APMISRR_totalTime_cost();
  // test_APMISRR_installment();
  // test_myAPMISRR_installment();
  // test_MISRRL_lambda();
  // test_MISRRL_installment();
  // test_MISRRL_all();
  // test_MISRRLL_lambda2();
  // test_MISRRLL_lambda1();
  // test_MISRRLL_installment();
  // test_MISRRLL_all();
  // test_MISRRLL_2_lambda();

  // compare_MISRR_and_MISRRL();

  return 0;
}