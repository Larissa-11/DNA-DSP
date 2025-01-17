# celery.py

# from __future__ import absolute_import, unicode_literals
# import os
# from celery import Celery
#
# # 设置默认的 Django settings 模块
# os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'djangoProject.settings')
#
# app = Celery('djangoProject')
#
# # 从 Django 的设置文件中加载配置
# app.config_from_object('django.conf:settings', namespace='CELERY')
#
# # 自动发现所有已安装的 Django app 中的 tasks.py 文件
# app.autodiscover_tasks()