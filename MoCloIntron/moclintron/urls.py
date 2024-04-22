from django.urls import path

from . import views

urlpatterns = [
    path("", views.index, name="index"),
    path("example", views.example, name="example"),
    path("results", views.results, name="results"),
    path("description", views.description, name="description"),
    path("double_slider", views.description, name="double_slider"),
    path("embed_test", views.embed_test, name="embed_test"),
    path("optimize_codons", views.optimize_codons, name="optimize_codons"),
    path("htmx_mouse_entered", views.htmx_mouse_entered, name="htmx_mouse_entered"),
    path("htmx_messages", views.htmx_messages, name="htmx_messages"),
    path("htmx_click", views.htmx_click, name="htmx_click"),
    path("htmx_click_site", views.htmx_click_site, name="htmx_click_site"),
    path("insert", views.insert, name="insert"),
       
]