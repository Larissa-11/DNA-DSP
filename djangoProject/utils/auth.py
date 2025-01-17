# -*- coding:utf-8 -*-
# @Time : 2023/11/25 21:12
# @Author: fbz
# @File : auth.py
from django.contrib.auth.mixins import LoginRequiredMixin
from django.http import JsonResponse


# 验证用户是否登陆
# LoginRequiredMixin未登錄是返回重定向，前後端分離需要返回json
class LoginRequiredJsonMixin(LoginRequiredMixin):
    pass
    # 重寫方法
    # def handle_no_permission(self):
    #     return JsonResponse({'code': 401, 'msg': '沒有登錄'}, status=401)
