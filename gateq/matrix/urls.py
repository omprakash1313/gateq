from . import views
from django.urls import path 

urlpatterns =[
    path('',views.index,name='home'),
    path('composer',views.composer,name='composer'),
    path('matAnimation',views.matAnimation,name='matAnimation'),
]