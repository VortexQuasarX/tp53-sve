import google.generativeai as genai
import sys

API_KEY = "AIzaSyDc7DEAAVxDFaSv3WQ5k0Knd3MmZtVweGQ"
genai.configure(api_key=API_KEY)

try:
    print("Available Models:")
    for m in genai.list_models():
        if 'generateContent' in m.supported_generation_methods:
            print(m.name)
except Exception as e:
    print(f"Error querying models: {e}")
