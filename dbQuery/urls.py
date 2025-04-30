from django.urls import path

from . import views

urlpatterns = [
    path("", views.home, name="home"),
    path('process_form/', views.process_form, name='process_form'),
    path("documentation/", views.documentation, name = "documentation"),
    path("dataset/", views.dataset, name = "dataset"),
    path("contact/", views.contact, name = "contact"),
    path("download_proteins/", views.download_proteins, name='download_proteins'),
    path('download_excel/', views.download_excel, name='download_excel'),
    path('base/',views.base, name='base'),
]