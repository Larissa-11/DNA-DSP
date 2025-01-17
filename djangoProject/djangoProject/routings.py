from django.urls import re_path

from djangoProject import consumers

websocket_urlpatterns = [
    # 示例 url ： xxxxx/room/x1/
    re_path('ws/notify/', consumers.NotificationConsumer.as_asgi()),
]
