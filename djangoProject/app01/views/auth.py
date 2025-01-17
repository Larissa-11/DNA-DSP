import json
import re

from django import forms
from django.contrib.auth.models import User
from django.http import JsonResponse
from django.shortcuts import redirect
from django.views import View
from django.contrib.auth import authenticate, login, logout
from django.core.handlers.wsgi import WSGIRequest

from djangoProject.consumers import logger


class UserRegistrationForm(forms.ModelForm):
    class Meta:
        model = User
        fields = ['username', 'password']

    def clean_password(self):
        cleaned_data = super().clean()
        password = cleaned_data.get('password')

        if len(password) != 8:
            raise forms.ValidationError("The password must be 8 digits")

        return password

    def clean_username(self):
        username = self.cleaned_data.get('username')
        if not re.match(r'^[a-zA-Z0-9_.+-]+@[a-zA-Z0-9-]+\.[a-zA-Z0-9-.]+$', username):
            raise forms.ValidationError("The username does not comply with email rules")
        if User.objects.filter(username=username).exists():
            raise forms.ValidationError("The user already exists")
        return username

    def save(self, commit=True):
        user = super().save(commit=False)
        user.set_password(self.cleaned_data["password"])
        if commit:
            user.save()
        return user


class RsgisterView(View):

    def post(self, request):
        try:
            data = json.loads(request.body.decode())
        except:
            return JsonResponse({'code': 400, 'msg': "Parameter error"})

        form = UserRegistrationForm(data)
        if form.is_valid():
            try:
                user = form.save()  # 将数据保存到数据库中
            except Exception as exc:
                return JsonResponse({'code': 400, 'msg': str(exc)})
            # 状态保持 -- 登录用户的状态保持
            login(request, user)
            return JsonResponse({'code': 200, 'msg': 'login was successful'})
        else:
            errors = json.loads(form.errors.as_json())
            msg = list(errors.values())[0][0]["message"]
            return JsonResponse({'code': 400, 'msg': msg})


class LoginView(View):

    def post(self, request: WSGIRequest):
        data = json.loads(request.body.decode())
        username = data.get('username')
        password = data.get('password')
        remembered = data.get('remembered')

        if not all([username, password]):
            return JsonResponse({'code': 400, 'msg': 'Incomplete parameters'})

        # 如果用户不存在提示注册
        # 强烈不推荐这个验证，因为可能造成攻击者确定用户名，建议直接在下面写账号或密码错误
        user = User.objects.filter(username=username).exists()
        if not user:
            return JsonResponse({'code': 400, 'msg': 'User does not exist, please register'})

        # authenticate傳遞用戶名和密碼
        # 如果用戶名和密碼正確，返回User信息，不正確返回None
        user = authenticate(username=username, password=password)
        if user is None:
            return JsonResponse({'code': 400, 'msg': 'Password input error'})

        login(request, user)
        # 是否免登录
        if remembered:
            # 14天免登录
            # set_expiry(value), value是一个整数，将在value秒后过期
            #                           None 默认两周，可以在setting中设置
            #                            0   将在浏览器关闭时过期
            request.session.set_expiry(None)
        else:
            # 不記住密碼，瀏覽器關閉session過期
            request.session.set_expiry(0)

        # 返回响应
        response = JsonResponse({'code': 200, 'msg': 'Login succeeded'})
        # 爲了首頁顯示用戶信息展示
        response.set_cookie('username', user.username)

        return response


class LogoutView(View):

    def post(self, request):
        # 刪除session信息
        logout(request)

        # return redirect('/login/')
        response = JsonResponse({'code': 200, 'msg': 'Successfully logged out and logged in'})
        response.delete_cookie('username')
        return response

class ModifyPasswordView(View):

    def post(self, request):
            old_pass = request.POST.get('old_pass')
            new_pass = request.POST.get('new_pass')

            try:
                request.user.check_password(old_pass)
            except Exception as e:
                logger.error(e)
                return JsonResponse({'code': 200, 'msg': "Parameter error"})

            # 修改密码
            try:
                request.user.set_password(new_pass)
                request.user.save()
                logout(request)
                return JsonResponse({'code': 200, 'msg': "success"})
            except Exception as e:
                logger.error(e)
                return JsonResponse({'code': 200, 'msg': "Parameter error"})

