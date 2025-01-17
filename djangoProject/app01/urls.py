# -*- coding:utf-8 -*-
# @Time : 2023/11/25 20:41
# @Author: fbz
# @File : urls.py
from django.urls import path

from app01.views import auth, coding_scheme, sr, algorithm_evaluation, reconstruction_analysis, errorSimulator, \
    primer_design

urlpatterns = [
    path("login/", auth.LoginView.as_view()),
    path("logout/", auth.LogoutView.as_view()),
    path("regiest/", auth.RsgisterView.as_view()),
    path("encode/", coding_scheme.EncodeView.as_view()),
    path("Decode/", coding_scheme.DecodeView.as_view()),
    path("sr/", sr.ClustringView.as_view()),
    path("Encode_Evalution/", algorithm_evaluation.AlgorithmEvaluationView.as_view()),  # 没有使用
    path("test/", algorithm_evaluation.Test.as_view()),  # 异步的algorithm_evaluation接口
    path("tests/", algorithm_evaluation.GetData.as_view()),  # 初始化时获取图表数据的接口，此接口无数据时，展示默认数据
    path("goto/", algorithm_evaluation.GetTask.as_view()),  # 初始化时获取消息接口
    path("reconstruction_analysis/", reconstruction_analysis.ClustringAnalysis.as_view()),  # 初始化时获取消息接口
    path("ErrorSimulator/", errorSimulator.ErrorSimulator.as_view()),  # error 计算
    path("primer_design/", primer_design.primer_design.as_view()),  # error 计算
    path("GetTaskDataInfo/", algorithm_evaluation.GetTaskDataInfo.as_view()),  # 获取消息内容
    path("ModifyPasswor/", auth.ModifyPasswordView.as_view()),  # 修改密码
]
