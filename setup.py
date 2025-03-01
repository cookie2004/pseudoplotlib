from setuptools import setup, find_packages

setup(
    name='pseudoplotlib',  # Choose a name for your package (e.g., pseudoplotlib)
    version='0.1.0',      # Start with a version number (e.g., 0.1.0)
    packages=find_packages(), # Automatically find all packages in your directory
    install_requires=[
        # List your dependencies here, e.g.,
        'numpy>=1.18',
        'matplotlib>=3.2',
    ],
    author='KittenMittens',     # Replace with your name
    author_email='cookie2004@gmail.com', # Replace with your email
    description='Matplotlib-based module for plotting pseudo protein representations',
    long_description=open('README.md').read(), # Read your README for long description
    long_description_content_type='text/markdown', # If your README is Markdown
    url='https://github.com/cookie2004/pseudoplotlib', # If you have a GitHub repo
    classifiers=[
        # Optional - Classifiers help users find your project
        'Development Status :: 3 - Alpha',  # Or "4 - Beta" or "5 - Production/Stable"
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',  # Replace with your license
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
    python_requires='>=3.6', # Specify the minimum Python version
)