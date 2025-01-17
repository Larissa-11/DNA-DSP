"""
ASGI config for djangoProject project.

It exposes the ASGI callable as a module-level variable named ``application``.

For more information on this file, see
https://docs.djangoproject.com/en/4.2/howto/deployment/asgi/
"""

import os

from channels.routing import ProtocolTypeRouter, URLRouter
from django.core.asgi import get_asgi_application

from djangoProject import routings

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "djangoProject.settings")

application = ProtocolTypeRouter({
    "http": get_asgi_application(),  # 自动找 urls.py ， 找视图函数  --》 http
    "websocket": URLRouter(routings.websocket_urlpatterns),  # routing(urls)、 consumers(views)
})
