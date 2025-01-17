# myproject/celery.py
from __future__ import absolute_import, unicode_literals
import os
from celery import Celery

from djangoProject import settings

# 设置环境变量为你的Django项目的settings模块
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'myproject.settings')

app = Celery('myproject')

# 从Django的settings.py模块中加载配置
app.config_from_object('django.conf:settings', namespace='CELERY')

# 自动发现每个app目录下的tasks.py模块
app.autodiscover_tasks(lambda: settings.INSTALLED_APPS)
