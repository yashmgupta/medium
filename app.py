# app.py

import streamlit as st
from data import research_summaries
from collections import defaultdict# app.py

import streamlit as st
from data import research_summaries
from collections import defaultdict

# -------------------------
# Set Page Configuration
# -------------------------
st.set_page_config(
    page_title="📚 Bioinformatics Learning App",
    layout="wide",
    initial_sidebar_state="expanded",
)

# -------------------------
# Custom CSS for Styling (Optional)
# -------------------------
def local_css(file_name):
    try:
        with open(file_name) as f:
            st.markdown(f"<style>{f.read()}</style>", unsafe_allow_html=True)
    except FileNotFoundError:
        st.warning("CSS file not found. Skipping custom styling.")

# Uncomment if you have a CSS file
# local_css("assets/styles.css")

# -------------------------
# Initialize Session State
# -------------------------
if 'current_index' not in st.session_state:
    st.session_state.current_index = 0

if 'likes' not in st.session_state:
    st.session_state.likes = defaultdict(bool)

if 'comments' not in st.session_state:
    st.session_state.comments = defaultdict(list)

if 'selected_category' not in st.session_state:
    st.session_state.selected_category = "All"

# -------------------------
# Sidebar: Category Selection
# -------------------------
st.sidebar.title("📖 Categories")

# Extract unique categories
categories = sorted(list(set([item['category'] for item in research_summaries])))

# Add 'All' option
categories.insert(0, "All")

# Single select for category
selected_category = st.sidebar.selectbox(
    "Select a Category",
    options=categories,
    index=0,
    key='category_select'
)

st.session_state.selected_category = selected_category

# -------------------------
# Main Title
# -------------------------
st.title("📚 Bioinformatics Learning App")

# -------------------------
# Search Functionality
# -------------------------
search_query = st.text_input(
    "🔍 Search Topics",
    value="",
    placeholder="Enter keywords to search...",
    help="Search within titles, summaries, and tags."
)

# -------------------------
# Filter Summaries Based on Category and Search
# -------------------------
def filter_summaries(summaries, category, query):
    filtered = summaries.copy()
    if category != "All":
        filtered = [s for s in filtered if s['category'] == category]
    if query:
        query = query.lower()
        filtered = [
            s for s in filtered
            if query in s['title'].lower()
            or query in s['content'].lower()
            or any(query in tag.lower() for tag in s['tags'])
        ]
    return filtered

filtered_summaries = filter_summaries(research_summaries, st.session_state.selected_category, search_query)

# -------------------------
# Handle No Results
# -------------------------
if not filtered_summaries:
    st.warning(f"No results found for '{search_query}' in category '{selected_category}'.")
    st.stop()

# -------------------------
# Pagination: Current Summary
# -------------------------
current_summary = filtered_summaries[st.session_state.current_index]

# -------------------------
# Display Summary
# -------------------------
with st.container():
    st.subheader(current_summary['title'])
    st.markdown(current_summary['content'])
    st.markdown(f"**Category:** {current_summary['category']}")
    st.markdown(f"**Tags:** {', '.join(current_summary['tags'])}")

    # -------------------------
    # Action Buttons
    # -------------------------
    col1, col2, col3 = st.columns(3)

    with col1:
        like_button = st.button(
            "❤️ Like",
            key=f"like_{current_summary['id']}",
            on_click=lambda: toggle_like(current_summary['id']),
        )
        if st.session_state.likes[current_summary['id']]:
            st.markdown("**You liked this!**")

    with col2:
        comment_input = st.text_input(
            "💬 Add a comment...",
            key=f"comment_input_{current_summary['id']}",
            placeholder="Type your comment here..."
        )
        if st.button("Submit Comment", key=f"submit_comment_{current_summary['id']}"):
            if comment_input:
                st.session_state.comments[current_summary['id']].append(comment_input)
                st.success("Comment submitted!")
            else:
                st.warning("Please enter a comment before submitting.")

    with col3:
        share_button = st.button("🔗 Share", key=f"share_{current_summary['id']}")
        if share_button:
            st.info("Share functionality is not implemented yet.")

    # -------------------------
    # Display Comments
    # -------------------------
    if st.session_state.comments[current_summary['id']]:
        st.markdown("### 💬 Comments")
        for idx, comment in enumerate(st.session_state.comments[current_summary['id']], 1):
            st.write(f"{idx}. {comment}")

# -------------------------
# Navigation Buttons
# -------------------------
st.markdown("---")
nav_col1, nav_col2, nav_col3 = st.columns([1, 2, 1])

with nav_col1:
    if st.button("⬅️ Previous", key="prev_button"):
        if st.session_state.current_index > 0:
            st.session_state.current_index -= 1
        else:
            st.warning("You're at the first topic.")

with nav_col3:
    if st.button("Next ➡️", key="next_button"):
        if st.session_state.current_index < len(filtered_summaries) - 1:
            st.session_state.current_index += 1
        else:
            st.warning("You've reached the last topic.")

# -------------------------
# Pagination Indicator
# -------------------------
st.markdown(
    f"**Topic {st.session_state.current_index + 1} of {len(filtered_summaries)}**"
)

# -------------------------
# Toggle Like Function
# -------------------------
def toggle_like(summary_id):
    st.session_state.likes[summary_id] = not st.session_state.likes[summary_id]

import pandas as pd
import plotly.express as px

# Set page configuration
st.set_page_config(
    page_title="📚 Bioinformatics Research Reels",
    layout="wide",
    initial_sidebar_state="expanded",
)

# Custom CSS for better styling (optional)
def local_css(file_name):
    try:
        with open(file_name) as f:
            st.markdown(f"<style>{f.read()}</style>", unsafe_allow_html=True)
    except FileNotFoundError:
        pass

# Uncomment if you have a styles.css file
# local_css("styles.css")

# Initialize session state variables
if 'current_index' not in st.session_state:
    st.session_state.current_index = 0

if 'likes' not in st.session_state:
    st.session_state.likes = defaultdict(bool)

if 'comments' not in st.session_state:
    st.session_state.comments = defaultdict(list)

# Sidebar for category selection
st.sidebar.title("📖 Categories")

# Get unique categories
categories = sorted(list(set([item['category'] for item in research_summaries])))

# Multi-select for categories
selected_categories = st.sidebar.multiselect(
    "Select Categories",
    options=categories,
    default=categories,  # All selected by default
)

# Main title
st.title("📚 Bioinformatics Research Reels")

# Search bar
search_query = st.text_input(
    "🔍 Search topics...",
    value="",
    placeholder="Enter keywords to search...",
)

# Function to filter summaries based on category and search
def filter_summaries(summaries, categories, query):
    filtered = [summary for summary in summaries if summary['category'] in categories]
    if query:
        query = query.lower()
        filtered = [
            summary for summary in filtered
            if query in summary['title'].lower()
            or query in summary['summary'].lower()
            or any(query in tag.lower() for tag in summary['tags'])
            or any(query in author.lower() for author in summary['authors'])
        ]
    return filtered

# Apply filters
filtered_summaries = filter_summaries(research_summaries, selected_categories, search_query)

# Display data visualization (e.g., Publications per Year)
def display_publications_per_year(summaries):
    df = pd.DataFrame(summaries)
    if df.empty:
        st.write("No data available for visualization.")
        return
    pub_counts = df.groupby('year').size().reset_index(name='publications')
    fig = px.bar(pub_counts, x='year', y='publications',
                 title='Number of Publications per Year',
                 labels={'year': 'Year', 'publications': 'Publications'})
    st.plotly_chart(fig, use_container_width=True)

# Display data visualization
display_publications_per_year(filtered_summaries)

# Handle case with no results
if not filtered_summaries:
    st.warning(f"No results found for '{search_query}' in selected categories.")
    st.stop()

# Pagination: Determine the current summary
current_summary = filtered_summaries[st.session_state.current_index]

# Display the summary in an expandable section
with st.container():
    # Summary Header
    summary_header_col1, summary_header_col2 = st.columns([3, 1])
    with summary_header_col1:
        st.subheader(current_summary['title'])
        st.write(f"**Authors:** {', '.join(current_summary['authors'])} | **Year:** {current_summary['year']}")
    with summary_header_col2:
        # Display like count or other metadata if available
        pass  # Placeholder for future metadata

    # Summary Content
    st.write(current_summary['summary'])
    st.write(f"**Tags:** {', '.join(current_summary['tags'])}")

    # Action buttons
    action_col1, action_col2, action_col3 = st.columns(3)

    with action_col1:
        like_button = st.button(
            "❤️ Like",
            key=f"like_{current_summary['id']}",
            on_click=lambda: toggle_like(current_summary['id']),
        )
        if st.session_state.likes[current_summary['id']]:
            st.markdown("**You liked this!**")

    with action_col2:
        comment_input = st.text_input(
            "💬 Add a comment...",
            key=f"comment_input_{current_summary['id']}",
            placeholder="Type your comment here..."
        )
        if st.button("Submit Comment", key=f"submit_comment_{current_summary['id']}"):
            if comment_input:
                st.session_state.comments[current_summary['id']].append(comment_input)
                st.success("Comment submitted!")
            else:
                st.warning("Please enter a comment before submitting.")

    with action_col3:
        share_button = st.button("🔗 Share", key=f"share_{current_summary['id']}")
        if share_button:
            st.info("🔗 Share functionality is not implemented yet.")

    # Display comments if any
    if st.session_state.comments[current_summary['id']]:
        st.markdown("### 💬 Comments")
        for idx, comment in enumerate(st.session_state.comments[current_summary['id']], 1):
            st.write(f"{idx}. {comment}")

# Navigation buttons
st.markdown("---")
nav_col1, nav_col2, nav_col3 = st.columns([1, 2, 1])

with nav_col1:
    if st.button("⬅️ Previous", key="prev_button"):
        if st.session_state.current_index > 0:
            st.session_state.current_index -= 1
        else:
            st.warning("You're at the first summary.")

with nav_col3:
    if st.button("Next ➡️", key="next_button"):
        if st.session_state.current_index < len(filtered_summaries) - 1:
            st.session_state.current_index += 1
        else:
            st.warning("You've reached the last summary.")

# Pagination Indicator
st.markdown(
    f"**Summary {st.session_state.current_index + 1} of {len(filtered_summaries)}**"
)

# Function to toggle like
def toggle_like(summary_id):
    st.session_state.likes[summary_id] = not st.session_state.likes[summary_id]
