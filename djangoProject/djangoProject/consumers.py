import json
from channels.generic.websocket import AsyncWebsocketConsumer
import logging

logger = logging.getLogger(__name__)

class NotificationConsumer(AsyncWebsocketConsumer):
    async def connect(self):
        self.group_name = 'notifications'

        # 加入通知组
        await self.channel_layer.group_add(
            self.group_name,
            self.channel_name
        )

        await self.accept()

    async def disconnect(self, close_code):
        # 离开通知组
        await self.channel_layer.group_discard(
            self.group_name,
            self.channel_name
        )

    async def receive(self, text_data):
        if text_data:
            try:
                tmp_data = json.loads(text_data)
                message = tmp_data['message']
                id = tmp_data['id']
                dataType = tmp_data['dataType']
                data = tmp_data.get("data")

                # 发送消息到组
                await self.channel_layer.group_send(
                    self.group_name,
                    {
                        'type': 'send_notification',
                        'message': message,
                        'id': id,
                        'dataType': dataType,
                        'data': data,
                    }
                )
            except json.JSONDecodeError as e:
                logger.error(f"JSON decode error: {e}")
                await self.send(text_data=json.dumps({
                    'error': 'Invalid JSON format'
                }))
        else:
            logger.error("Received empty text_data")

    async def send_notification(self, event):
        message = event['message']
        id = event['id']
        dataType = event['dataType']
        data = event['data']
        error = event['error']

        # 发送消息到WebSocket
        await self.send(text_data=json.dumps({
            'message': message,
            'id': id,
            'error': error,
            'dataType': dataType,
            'data': data,
        }))
